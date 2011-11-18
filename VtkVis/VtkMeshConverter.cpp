/**
 * \file VtkMeshConverter.cpp
 * 23/08/2011 KR Initial implementation
 *
 */

#include "VtkMeshConverter.h"

#include "GridAdapter.h"

// Conversion from Image to QuadMesh
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

// Conversion from vtkUnstructuredGrid
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>

MeshLib::CFEMesh* VtkMeshConverter::convertImgToMesh(vtkImageData* img,
                                                     const std::pair<double,double> &origin,
                                                     const double &scalingFactor,
													 MshElemType::type t,
													 bool setAsElevation)
{
	if ((t != MshElemType::TRIANGLE) && (t != MshElemType::QUAD))
	{
		std::cout << "Error in VtkMeshConverter::convertImgToMesh() - Invalid Mesh Element Type..." << std::endl;
		return NULL;
	}

	vtkSmartPointer<vtkUnsignedCharArray> pixelData = vtkSmartPointer<vtkUnsignedCharArray>(
	        vtkUnsignedCharArray::SafeDownCast(img->GetPointData()->GetScalars()));
	int* dims = img->GetDimensions();

	MeshLib::CFEMesh* mesh(new MeshLib::CFEMesh());
	const size_t imgHeight = dims[0];
	const size_t imgWidth  = dims[1];
	double* pixVal (new double[imgHeight * imgWidth]);
	bool* visNodes(new bool[imgWidth * imgHeight]);

	for (size_t i = 0; i < imgWidth; i++)
		for (size_t j = 0; j < imgHeight; j++)
		{
			const size_t index = i * imgHeight + j;
			const double* colour = pixelData->GetTuple4(index);
			pixVal[index] = 0.3 * colour[0] + 0.6 * colour[1] + 0.1 * colour[2];
			double zValue = (setAsElevation) ? pixVal[index] : 0.0;
			double coords[3] = { origin.first + (scalingFactor * j), 
								 origin.second + (scalingFactor * i), 
								 zValue };
			visNodes[index] = (colour[3] > 0);

			MeshLib::CNode* node(new MeshLib::CNode(index));
			node->SetCoordinates(coords);
			mesh->nod_vector.push_back(node);
		}

	// set mesh elements
	for (size_t i = 0; i < imgWidth - 1; i++)
		for (size_t j = 0; j < imgHeight - 1; j++)
		{
			const int index = i * imgHeight + j;

			// if node is visible
			if (visNodes[index])
			{
				const int mat = (setAsElevation) ? 0 : static_cast<int>(pixVal[index]);
				if (t == MshElemType::TRIANGLE)
				{
					mesh->ele_vector.push_back(createElement(t, mat, index, index + 1, 
															 index + imgHeight));       // upper left triangle
					mesh->ele_vector.push_back(createElement(t, mat, index + 1, 
															 index + imgHeight + 1, 
															 index + imgHeight));                   // lower right triangle
				}
				if (t == MshElemType::QUAD)
				{
					mesh->ele_vector.push_back(createElement(t, mat, index, index + 1, 
															 index + imgHeight + 1,
															 index + imgHeight));
				}
			}
		}

	mesh->ConstructGrid();
	delete [] pixVal;
	delete [] visNodes;
	return mesh;
}

MeshLib::CElem* VtkMeshConverter::createElement(MshElemType::type t, int mat, size_t node1, size_t node2, size_t node3, size_t node4)
{
	MeshLib::CElem* elem(new MeshLib::CElem());
	const size_t nNodes = (t == MshElemType::QUAD) ? 4 : 3;
	elem->setElementProperties(t);
	elem->SetPatchIndex(mat);
	elem->SetNodesNumber(nNodes);
	elem->SetNodeIndex(0, node1);
	elem->SetNodeIndex(1, node2);
	elem->SetNodeIndex(2, node3);
	if (t ==  MshElemType::QUAD)
		elem->SetNodeIndex(3, node4);
	elem->InitializeMembers();
	return elem;
}

MeshLib::CFEMesh* VtkMeshConverter::convertUnstructuredGrid(vtkUnstructuredGrid* grid)
{
	if (!grid)
		return NULL;

	MeshLib::CFEMesh* mesh(new MeshLib::CFEMesh());

	const size_t nNodes = grid->GetPoints()->GetNumberOfPoints();
	const size_t nElems = grid->GetNumberOfCells();

	// set mesh nodes
	double* coords = NULL;
	for (size_t i = 0; i < nNodes; i++)
	{
		coords = grid->GetPoints()->GetPoint(i);
		MeshLib::CNode* node(new MeshLib::CNode(i, coords[0], coords[1], coords[2]));
		mesh->nod_vector.push_back(node);
	}

	// set mesh elements
	vtkCell* cell(NULL);
	vtkDataArray* scalars = grid->GetCellData()->GetScalars("MaterialIDs");
	for (size_t i = 0; i < nElems; i++)
	{
		MeshLib::CElem* elem(new MeshLib::CElem());

		MshElemType::type elem_type = MshElemType::INVALID;
		int cell_type = grid->GetCellType(i);

		switch (cell_type)
		{
		case VTK_TRIANGLE:      elem_type = MshElemType::TRIANGLE;
			break;
		case VTK_QUAD:          elem_type = MshElemType::QUAD;
			break;
		case VTK_TETRA:         elem_type = MshElemType::TETRAHEDRON;
			break;
		case VTK_HEXAHEDRON:    elem_type = MshElemType::HEXAHEDRON;
			break;
		case VTK_WEDGE:         elem_type = MshElemType::PRISM;
			break;
		}

		if (elem_type != MshElemType::INVALID)
		{
			//elem->SetElementType(elem_type);
			elem->setElementProperties(elem_type);
			if (scalars)
				elem->SetPatchIndex(static_cast<int>(scalars->GetComponent(i,0)));  // HACK the name of the correct scalar array of the vtk file should probably be passed as an argument?!
		}
		else
		{
			std::cout <<
			"Error in GridAdapter::convertUnstructuredGrid() - Unknown mesh element type ..."
			          << std::endl;
			return NULL;
		}

		cell = grid->GetCell(i);
		size_t nElemNodes = cell->GetNumberOfPoints();
		elem->SetNodesNumber(nElemNodes);
		elem->nodes_index.resize(nElemNodes);

		for (size_t j = 0; j < nElemNodes; j++)
			elem->SetNodeIndex(j, cell->GetPointId(j));

		mesh->ele_vector.push_back(elem);
	}
	mesh->ConstructGrid();
	return mesh;
}
