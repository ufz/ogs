/**
 * \file GridAdapter.cpp
 * 24/03/2010 KR Initial implementation
 *
 */


#include "GridAdapter.h"

#include <iostream>
#include <list>
#include <cstdlib>
#include "StringTools.h"

// Conversion from Image to QuadMesh
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

// Conversion from vtkUnstructuredGrid
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkCellData.h>

using Mesh_Group::CFEMesh;

GridAdapter::GridAdapter(const Mesh_Group::CFEMesh* mesh) :
	_name(""), _nodes(new std::vector<GEOLIB::Point*>), _elems(new std::vector<Element*>), _mesh(mesh)
{
	if (mesh) convertCFEMesh(mesh);
}

GridAdapter::GridAdapter(const std::string &filename) :
	_name(""), _nodes(new std::vector<GEOLIB::Point*>), _elems(new std::vector<Element*>), _mesh(NULL)
{
	readMeshFromFile(filename);
}

GridAdapter::~GridAdapter()
{
	size_t nNodes = _nodes->size();
	size_t nElems = _elems->size();

	for (size_t i=0; i<nNodes; i++) delete (*_nodes)[i];
	for (size_t i=0; i<nElems; i++) delete (*_elems)[i];

	delete this->_nodes;
	delete this->_elems;
}


int GridAdapter::convertCFEMesh(const Mesh_Group::CFEMesh* mesh)
{
	if (!mesh) return 0;
	
	size_t nNodes = mesh->nod_vector.size();
	for (size_t i=0; i<nNodes; i++)
	{
		GEOLIB::Point* pnt = new GEOLIB::Point(mesh->nod_vector[i]->X(), mesh->nod_vector[i]->Y(), mesh->nod_vector[i]->Z());
		_nodes->push_back(pnt);
	}

	Element* newElem = NULL;
	size_t nElems = mesh->ele_vector.size();
	size_t nElemNodes = 0;

	for (size_t i=0; i<nElems; i++)
	{
		newElem = new Element();
		newElem->type = mesh->ele_vector[i]->GetElementType();
		newElem->material = mesh->ele_vector[i]->GetPatchIndex();

		if (newElem->type != MshElemType::INVALID)
		{
			std::vector<long> elemNodes;
			nElemNodes = mesh->ele_vector[i]->nodes_index.Size();
			for (size_t j=0; j<nElemNodes; j++)
				newElem->nodes.push_back(mesh->ele_vector[i]->GetNodeIndex(j));

			_elems->push_back(newElem);
		}
		else
			std::cout << "GridAdapter::convertCFEMesh() - Error recognising element type..." << std::endl;
	}

	return 1;
}

const CFEMesh* GridAdapter::getCFEMesh() const
{
	if (_mesh) return _mesh;
	return toCFEMesh();
}

const CFEMesh* GridAdapter::toCFEMesh() const
{
	CFEMesh* mesh (new CFEMesh());
	std::cout << "Converting mesh object ... ";

	// set mesh nodes
	size_t nNodes = _nodes->size();
	for (size_t i=0; i<nNodes; i++)
	{
		Mesh_Group::CNode* node(new Mesh_Group::CNode(i));
		double coords[3] = { (*(*_nodes)[i])[0], (*(*_nodes)[i])[1], (*(*_nodes)[i])[2] };
		node->SetCoordinates(coords);
		mesh->nod_vector.push_back(node);
	}

	// set mesh elements
	size_t nElems = _elems->size();
	for (size_t i=0; i<nElems; i++)
	{
		Mesh_Group::CElem* elem(new Mesh_Group::CElem());
		elem->SetElementType((*_elems)[i]->type);
		elem->SetPatchIndex((*_elems)[i]->material);

		size_t nElemNodes = ((*_elems)[i]->nodes).size();
		elem->SetNodesNumber(nElemNodes);
		elem->nodes_index.resize(nElemNodes);
		for (size_t j=0; j<nElemNodes; j++)
		{
			int a = (*_elems)[i]->nodes[j];
			//elem->nodes_index[j] = a;
			elem->SetNodeIndex(j, a);
		}

		mesh->ele_vector.push_back(elem);
	}
	// nfaces?
	// nedges?
	mesh->ConstructGrid();
	std::cout << "done." << std::endl;

	return mesh;
}

int GridAdapter::readMeshFromFile(const std::string &filename)
{
	std::string line;
	std::list<std::string>::const_iterator it;

	std::ifstream in( filename.c_str() );

	if (!in.is_open())
	{
		std::cout << "GridAdapter::readMeshFromFile() - Could not open file..."  << std::endl;
		return 0;
	}

	// try to find the start of the nodes list
	while ( getline(in, line) )
	{
		trim(line);
		if (line.compare("$NODES") == 0) break;
	}

	// read number of nodes
	getline(in, line);
	trim(line);

	// read all nodes
	while ( getline(in, line) )
	{
		trim(line);
		if (line.compare("$ELEMENTS") == 0) break;

		std::list<std::string> fields = splitString(line, ' ');

		if (fields.size() >= 4)
		{
			it = fields.begin();

			if (atoi(it->c_str()) == (int)_nodes->size())
			{
				GEOLIB::Point* pnt = new GEOLIB::Point();

				(*pnt)[0] = strtod((++it)->c_str(), 0);
				(*pnt)[1] = strtod((++it)->c_str(), 0);
				(*pnt)[2] = strtod((++it)->c_str(), 0);

				_nodes->push_back(pnt);
			}
			else
				std::cout << "GridAdapter::readMeshFromFile() - Index error while reading nodes..." << std::endl;
		}
		else
			std::cout << "GridAdapter::readMeshFromFile() - Error reading node format..." << std::endl;
	}

	if (line.compare("$ELEMENTS") == 0)
	{
		Element* newElem;

		// read number of elements
		getline(in, line);
		trim(line);

		// read all elements
		while ( getline(in, line) )
		{
			trim(line);
			if (line.compare("$LAYER") == 0) break;

			std::list<std::string> fields = splitString(line, ' ');

			if (fields.size() >= 6) {

				it = fields.begin();
				if (atoi(it->c_str()) == (int)_elems->size())
				{
					newElem = new Element();

					if ((++it)->empty()) it++;
					newElem->material = atoi(it->c_str());	// material group
					if ((++it)->empty()) it++;
					newElem->type = getElementType(*it);	// element type

					if (newElem->type != MshElemType::INVALID)
					{
						while ((++it) != fields.end())
						{
							if (it->empty()) continue;
							newElem->nodes.push_back(atoi(it->c_str()));	// next node id
						}

						_elems->push_back(newElem);
					}
					else
						std::cout << "GridAdapter::readMeshFromFile() - Error recognising element type..." << std::endl;
				}
				else
					std::cout << "GridAdapter::readMeshFromFile() - Index error while reading elements..." << std::endl;
			}
			else
				std::cout << "GridAdapter::readMeshFromFile() - Error reading element format..." << std::endl;
		}
	}
	else
		std::cout << "GridAdapter::readMeshFromFile() - Index error after reading nodes..." << std::endl;

	in.close();

	return 1;
}

MshElemType::type GridAdapter::getElementType(const std::string &t) const
{
	if (t.compare("tri") == 0)  return MshElemType::TRIANGLE;
	if (t.compare("line") == 0) return MshElemType::LINE;
	if (t.compare("quad") == 0) return MshElemType::QUAD;
	if (t.compare("tet") == 0)  return MshElemType::TETRAHEDRON;
	if (t.compare("hex") == 0)  return MshElemType::HEXAHEDRON;
	if (t.compare("pri") == 0)  return MshElemType::PRISM;
	else return MshElemType::INVALID;
}

size_t GridAdapter::getNumberOfMaterials() const
{
	size_t nElems = _elems->size();
	size_t maxMaterialID = 0;

	for (size_t i=0; i<nElems; i++)
	{
		if ((*_elems)[i]->material > maxMaterialID) maxMaterialID = (*_elems)[i]->material;
	}
	return maxMaterialID;
}

const std::vector<GridAdapter::Element*> *GridAdapter::getElements(size_t matID) const
{
	std::vector<GridAdapter::Element*> *matElems = new std::vector<GridAdapter::Element*>;
	size_t nElements = _elems->size();
	for (size_t i=0; i<nElements; i++)
	{
		if ((*_elems)[i]->material == matID)
			matElems->push_back((*_elems)[i]);
	}

	return matElems;
}

Mesh_Group::CFEMesh* GridAdapter::convertImgToMesh(vtkImageData* img, const std::pair<double,double> &origin, const double &scalingFactor)
{
	vtkSmartPointer<vtkUnsignedCharArray> pixelData = vtkSmartPointer<vtkUnsignedCharArray>(vtkUnsignedCharArray::SafeDownCast(img->GetPointData()->GetScalars()));
	int* dims = img->GetDimensions();

	Mesh_Group::CFEMesh* mesh(new CFEMesh());
	size_t imgHeight = dims[0];
	size_t imgWidth  = dims[1];
	std::vector<size_t> visNodes(imgWidth*imgHeight);

	for (size_t i=0; i<imgWidth; i++)
	{
		for (size_t j=0; j<imgHeight; j++)
		{
			size_t index = i*imgHeight+j;
			const double* colour = pixelData->GetTuple4(index);
			double pixcol = 0.3*colour[0] + 0.6*colour[1] + 0.1*colour[2];
			double coords[3] = { origin.first+(scalingFactor*j), origin.second+(scalingFactor*i), pixcol };
			visNodes[index] = (colour[3]>0);

			Mesh_Group::CNode* node(new Mesh_Group::CNode(index));
			node->SetCoordinates(coords);
			mesh->nod_vector.push_back(node);
		}
	}

	// set mesh elements
	for (size_t i=0; i<imgWidth-1; i++)
    {
		for (size_t j=0; j<imgHeight-1; j++)
		{
			int index = i*imgHeight+j;

			// if node is visible
			if (visNodes[index])
			{

				mesh->ele_vector.push_back(createElement(index, index+1, index+imgHeight)); // upper left triangle
				mesh->ele_vector.push_back(createElement(index+1, index+imgHeight+1, index+imgHeight)); // lower right triangle
			}
		}
	}
	mesh->ConstructGrid();
	return mesh;
}

Mesh_Group::CElem* GridAdapter::createElement(size_t node1, size_t node2, size_t node3)
{
	Mesh_Group::CElem* elem(new Mesh_Group::CElem());
	elem->setElementProperties(MshElemType::TRIANGLE);
	elem->SetPatchIndex(1);
	elem->SetNodesNumber(3);
	elem->SetNodeIndex(0, node1);
	elem->SetNodeIndex(1, node2);
	elem->SetNodeIndex(2, node3);
	elem->InitializeMembers();
	return elem;
}

Mesh_Group::CFEMesh* GridAdapter::convertUnstructuredGrid(vtkUnstructuredGrid* grid)
{
	if (!grid) return NULL;

	Mesh_Group::CFEMesh* mesh(new Mesh_Group::CFEMesh());

	size_t nNodes = grid->GetPoints()->GetNumberOfPoints();
	size_t nElems = grid->GetNumberOfCells();

	// set mesh nodes
	double* coords = NULL;
	for (size_t i=0; i<nNodes; i++)
	{
		coords = grid->GetPoints()->GetPoint(i);
		Mesh_Group::CNode* node(new Mesh_Group::CNode(i, coords[0], coords[1], coords[2]));
		mesh->nod_vector.push_back(node);
	}
	
	// set mesh elements
	vtkCell* cell(NULL);
	vtkDataArray* scalars = grid->GetCellData()->GetScalars("MaterialIDs");
	for (size_t i=0; i<nElems; i++)
	{
		Mesh_Group::CElem* elem(new Mesh_Group::CElem());

		MshElemType::type elem_type = MshElemType::INVALID;
		int cell_type = grid->GetCellType(i);

		switch (cell_type)
		{
			case VTK_TRIANGLE:		elem_type=MshElemType::TRIANGLE;	break;
			case VTK_QUAD:			elem_type=MshElemType::QUAD;		break;
			case VTK_TETRA:			elem_type=MshElemType::TETRAHEDRON;	break;
			case VTK_HEXAHEDRON:	elem_type=MshElemType::HEXAHEDRON;	break;
			case VTK_WEDGE:			elem_type=MshElemType::PRISM;		break;
		}

		if (elem_type != MshElemType::INVALID)
			elem->SetElementType(elem_type);
			if (scalars) elem->SetPatchIndex(static_cast<int>(scalars->GetComponent(i,0))); // HACK the name of the correct scalar array of the vtk file should probably be passed as an argument?!
		else
		{
			std::cout << "Error in GridAdapter::convertUnstructuredGrid() - Unknown mesh element type ..." << std::endl;
			return NULL;
		}

		cell = grid->GetCell(i);
		size_t nElemNodes = cell->GetNumberOfPoints();
		elem->SetNodesNumber(nElemNodes);
		elem->nodes_index.resize(nElemNodes);

		for (size_t j=0; j<nElemNodes; j++)
		{
			elem->SetNodeIndex(j, cell->GetPointId(j));
		}

		mesh->ele_vector.push_back(elem);
	}
	mesh->ConstructGrid();
	return mesh;
}
