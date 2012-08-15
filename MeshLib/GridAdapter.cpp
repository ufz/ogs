/**
 * \file GridAdapter.cpp
 * 24/03/2010 KR Initial implementation
 *
 */

#include "GridAdapter.h"

#include "StringTools.h"
#include <cstdlib>
#include <iostream>
#include <list>

GridAdapter::GridAdapter(const MeshLib::CFEMesh* mesh) :
	_name(""), _nodes(new std::vector<GEOLIB::Point*>), _elems(new std::vector<Element*>),
	_mesh(mesh)
{
	if (mesh)
		this->convertCFEMesh(mesh);
}

GridAdapter::GridAdapter(const std::string &filename) :
	_name(""), _nodes(new std::vector<GEOLIB::Point*>), _elems(new std::vector<Element*>),
	_mesh(NULL)
{
	readMeshFromFile(filename);
}

GridAdapter::GridAdapter(const GridAdapter* grid) :
	_name(""), _nodes(new std::vector<GEOLIB::Point*>), _elems(new std::vector<Element*>),
	_mesh(NULL)
{
	if (grid)
	{
		this->_name = grid->getName();
	
		const std::vector<GEOLIB::Point*> *nodes = grid->getNodes();
		const size_t nNodes(nodes->size());
		for (size_t i=0; i<nNodes; i++)
		{
			GEOLIB::Point* pnt (new GEOLIB::Point((*nodes)[i]->getData()));
			this->addNode(pnt);
		}

		const std::vector<GridAdapter::Element*> *elements = grid->getElements();
		const size_t nElems(elements->size());
		for (size_t i=0; i<nElems; i++)
		{
			GridAdapter::Element* elem = new GridAdapter::Element;
			elem->material = (*elements)[i]->material;
			const size_t nElemNodes ((*elements)[i]->nodes.size());
			for (size_t j=0; j<nElemNodes; j++)
				elem->nodes.push_back((*elements)[i]->nodes[j]);
			elem->type = (*elements)[i]->type;
			this->addElement(elem);
		}
	}	
}

GridAdapter::~GridAdapter()
{
	size_t nNodes = _nodes->size();
	size_t nElems = _elems->size();

	for (size_t i = 0; i < nNodes; i++)
		delete (*_nodes)[i];
	for (size_t i = 0; i < nElems; i++)
		delete (*_elems)[i];

	delete this->_nodes;
	delete this->_elems;
}

int GridAdapter::convertCFEMesh(const MeshLib::CFEMesh* mesh)
{
	if (!mesh)
		return 0;

	size_t nNodes = mesh->nod_vector.size();
	for (size_t i = 0; i < nNodes; i++)
		_nodes->push_back(new GEOLIB::Point(mesh->nod_vector[i]->getData()));

	size_t nElems = mesh->ele_vector.size();
	size_t nElemNodes = 0;

	for (size_t i = 0; i < nElems; i++)
	{
		Element* newElem = new Element();
		newElem->type = mesh->ele_vector[i]->GetElementType();
		newElem->material = mesh->ele_vector[i]->GetPatchIndex();

		if (newElem->type != MshElemType::INVALID)
		{
			std::vector<long> elemNodes;
			nElemNodes = mesh->ele_vector[i]->getNodeIndices().Size();
			for (size_t j = 0; j < nElemNodes; j++)
				newElem->nodes.push_back(mesh->ele_vector[i]->GetNodeIndex(j));

			_elems->push_back(newElem);
		}
		else
			std::cout << "GridAdapter::convertCFEMesh() - Error recognising element type..." << std::endl;
	}

	return 1;
}

const MeshLib::CFEMesh* GridAdapter::getCFEMesh() const
{
	if (_mesh)
		return _mesh;
	return toCFEMesh();
}

const MeshLib::CFEMesh* GridAdapter::toCFEMesh() const
{
	MeshLib::CFEMesh* mesh (new MeshLib::CFEMesh());
	std::cout << "Converting mesh object ... ";

	// set mesh nodes
	size_t nNodes = _nodes->size();
	for (size_t i = 0; i < nNodes; i++)
		mesh->nod_vector.push_back(new MeshLib::CNode(i, (*(*_nodes)[i])[0], (*(*_nodes)[i])[1], (*(*_nodes)[i])[2]));

	// set mesh elements
	size_t nElems = _elems->size();
	for (size_t i = 0; i < nElems; i++)
	{
		MeshLib::CElem* elem(new MeshLib::CElem());
		elem->setElementProperties((*_elems)[i]->type);
		elem->SetPatchIndex((*_elems)[i]->material);

		size_t nElemNodes = ((*_elems)[i]->nodes).size();
		for (size_t j = 0; j < nElemNodes; j++)
			elem->SetNodeIndex(j, (*_elems)[i]->nodes[j]);

		mesh->ele_vector.push_back(elem);
	}

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
		std::cout << "GridAdapter::readMeshFromFile() - Could not open file..."  <<
		std::endl;
		return 0;
	}

	// try to find the start of the nodes list
	while ( getline(in, line) )
	{
		trim(line);
		if (line.compare("$NODES") == 0)
			break;
	}

	// read number of nodes
	getline(in, line);
	trim(line);

	// read all nodes
	while ( getline(in, line) )
	{
		trim(line);
		if (line.compare("$ELEMENTS") == 0)
			break;

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
				std::cout <<
				"GridAdapter::readMeshFromFile() - Index error while reading nodes..."
				          << std::endl;
		}
		else
			std::cout <<
			"GridAdapter::readMeshFromFile() - Error reading node format..." <<
			std::endl;
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
			if (line.compare("$LAYER") == 0)
				break;

			std::list<std::string> fields = splitString(line, ' ');

			if (fields.size() >= 5)
			{
				it = fields.begin();
				if (atoi(it->c_str()) == (int)_elems->size())
				{
					newElem = new Element();

					if ((++it)->empty())
						it++;
					newElem->material = atoi(it->c_str()); // material group
					if ((++it)->empty())
						it++;
					newElem->type = getElementType(*it); // element type

					if (newElem->type != MshElemType::INVALID)
					{
						while ((++it) != fields.end())
						{
							if (it->empty())
								continue;
							newElem->nodes.push_back(atoi(it->c_str())); // next node id
						}

						_elems->push_back(newElem);
					}
					else
						std::cout <<
						"GridAdapter::readMeshFromFile() - Error recognising element type..."
						          << std::endl;
				}
				else
					std::cout <<
					"GridAdapter::readMeshFromFile() - Index error while reading element " << *it << "... "
					          << std::endl;
			}
			else
				std::cout <<
				"GridAdapter::readMeshFromFile() - Error reading element format..."
				          << std::endl;
		}
	}
	else
		std::cout <<
		"GridAdapter::readMeshFromFile() - Index error after reading nodes..." << std::endl;

	in.close();

	return 1;
}

MshElemType::type GridAdapter::getElementType(const std::string &t) const
{
	if (t.compare("tri") == 0)
		return MshElemType::TRIANGLE;
	if (t.compare("line") == 0)
		return MshElemType::LINE;
	if (t.compare("quad") == 0)
		return MshElemType::QUAD;
	if (t.compare("tet") == 0)
		return MshElemType::TETRAHEDRON;
	if (t.compare("hex") == 0)
		return MshElemType::HEXAHEDRON;
	if ((t.compare("pri") == 0) || (t.compare("pris") == 0))
		return MshElemType::PRISM;
	else
		return MshElemType::INVALID;
}

size_t GridAdapter::getNumberOfMaterials() const
{
	size_t nElems = _elems->size();
	size_t maxMaterialID = 0;

	for (size_t i = 0; i < nElems; i++)
		if ((*_elems)[i]->material > maxMaterialID)
			maxMaterialID = (*_elems)[i]->material;
	return maxMaterialID;
}

const std::vector<GridAdapter::Element*>* GridAdapter::getElements(size_t matID) const
{
	std::vector<GridAdapter::Element*>* matElems = new std::vector<GridAdapter::Element*>;
	size_t nElements = _elems->size();
	for (size_t i = 0; i < nElements; i++)
		if ((*_elems)[i]->material == matID)
			matElems->push_back((*_elems)[i]);

	return matElems;
}
