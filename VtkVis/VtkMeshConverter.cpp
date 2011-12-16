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
													 MshElemType::type elem_type,
													 UseIntensityAs::type intensity_type)
{
	if ((elem_type != MshElemType::TRIANGLE) && (elem_type != MshElemType::QUAD))
	{
		std::cout << "Error in VtkMeshConverter::convertImgToMesh() - Invalid Mesh Element Type..." << std::endl;
		return NULL;
	}

	vtkSmartPointer<vtkUnsignedCharArray> pixelData = vtkSmartPointer<vtkUnsignedCharArray>(
	        vtkUnsignedCharArray::SafeDownCast(img->GetPointData()->GetScalars()));
	int* dims = img->GetDimensions();

	const size_t imgHeight = dims[0];
	const size_t imgWidth  = dims[1];
	const size_t incHeight = imgHeight+1;
	const size_t incWidth  = imgWidth+1;
	double* pixVal (new double[incHeight * incWidth]);
	bool* visNodes(new bool[incWidth * incHeight]);
	int* node_idx_map(new int[incWidth * incHeight]);

	for (size_t j = 0; j < incWidth; j++)
	{
		pixVal[j]=0;
		visNodes[j]=false;
		node_idx_map[j]=-1;
	}
	for (size_t i = 0; i < imgWidth; i++)
	{
		for (size_t j = 0; j < imgHeight; j++)
		{
			const size_t img_idx = i * imgHeight + j;
			const size_t index = (i+1) * incHeight + j;
			const double* colour = pixelData->GetTuple4(img_idx);
			pixVal[index] = 0.3 * colour[0] + 0.6 * colour[1] + 0.1 * colour[2];
			visNodes[index] = (colour[3] > 0);
			node_idx_map[index]=-1;
		}
		pixVal[(i+2)*incHeight-1]=0;
		visNodes[(i+2)*incHeight-1]=false;
		node_idx_map[(i+2)*incHeight-1]=-1;
	}
	
	MeshLib::CFEMesh* mesh = constructMesh(pixVal, node_idx_map, visNodes, origin, imgHeight, imgWidth, scalingFactor, elem_type, intensity_type);

	delete [] pixVal;
	delete [] visNodes;
	delete [] node_idx_map;

	return mesh;
}

MeshLib::CFEMesh* VtkMeshConverter::convertImgToMesh(const double* img,
													 const std::pair<double,double> &origin,
													 const double imgHeight,
													 const double imgWidth,
													 const double &scalingFactor,
													 MshElemType::type elem_type,
													UseIntensityAs::type intensity_type)
{
	const size_t incHeight = imgHeight+1;
	const size_t incWidth  = imgWidth+1;
	double* pixVal (new double[incHeight * incWidth]);
	bool* visNodes(new bool[incWidth * incHeight]);
	int* node_idx_map(new int[incWidth * incHeight]);

	for (size_t j = 0; j < incWidth; j++)
	{
		pixVal[j]=0;
		visNodes[j]=false;
		node_idx_map[j]=-1;
	}
	for (size_t i = 0; i < imgWidth; i++)
	{
		for (size_t j = 0; j < imgHeight; j++)
		{
			const size_t img_idx = i * imgHeight + j;
			const size_t index = (i+1) * incHeight + j;
			pixVal[index] = img[img_idx];
			visNodes[index] = 1;
			node_idx_map[index]=-1;
		}
		pixVal[(i+2)*incHeight-1]=0;
		visNodes[(i+2)*incHeight-1]=false;
		node_idx_map[(i+2)*incHeight-1]=-1;
	}
	
	MeshLib::CFEMesh* mesh = constructMesh(pixVal, node_idx_map, visNodes, origin, imgHeight, imgWidth, scalingFactor, elem_type, intensity_type);

	delete [] pixVal;
	delete [] visNodes;
	delete [] node_idx_map;

	return mesh;
}

MeshLib::CFEMesh* VtkMeshConverter::constructMesh(const double* pixVal,
												  int* node_idx_map,
												  const bool* visNodes,
												  const std::pair<double,double> &origin,
                                                  const size_t &imgHeight,
												  const size_t &imgWidth,
                                                  const double &scalingFactor,
										 		  MshElemType::type elem_type,
												  UseIntensityAs::type intensity_type)
{
	const size_t incHeight = imgHeight+1;
	const size_t incWidth  = imgWidth+1;
	MeshLib::CFEMesh* mesh(new MeshLib::CFEMesh());
	size_t node_idx_count(0);
	const double x_offset(origin.first - scalingFactor/2.0);
	const double y_offset(origin.second - scalingFactor/2.0);

	for (size_t i = 0; i < incWidth; i++)
		for (size_t j = 0; j < incHeight; j++)
		{
			const size_t index = i * incHeight + j;

			bool set_node (false);
			if (j==0 && i==imgWidth) set_node = visNodes[index];
			else if (j==0)			 set_node = (visNodes[index] || visNodes[index+incHeight]);
			else if (i==imgWidth)	 set_node = (visNodes[index] || visNodes[index-1]);
			else					 set_node = (visNodes[index] || visNodes[index-1] || visNodes[index+incHeight] || visNodes[index+incHeight-1]);

			if (set_node)
			{
				double zValue = (intensity_type == UseIntensityAs::ELEVATION) ? pixVal[index] : 0.0;
				const double coords[3] = { x_offset + (scalingFactor * j),
									       y_offset + (scalingFactor * i),
									       zValue };

				MeshLib::CNode* node(new MeshLib::CNode(node_idx_count));
				node->SetCoordinates(coords);
				mesh->nod_vector.push_back(node);
				node_idx_map[index] = node_idx_count;
				node_idx_count++;
			}
		}

	// set mesh elements
	for (size_t i = 0; i < imgWidth; i++)
		for (size_t j = 0; j < imgHeight; j++)
		{
			const int index = i * incHeight + j;
			if ((node_idx_map[index]!=-1) && (node_idx_map[index+1]!=-1) && (node_idx_map[index+incHeight]!=-1) && (node_idx_map[index+incHeight+1]!=-1) && (visNodes[index+incHeight]))
			{
				const int mat = (intensity_type != UseIntensityAs::MATERIAL) ? 0 : static_cast<int>(pixVal[index+incHeight]);
				if (elem_type == MshElemType::TRIANGLE)
				{
					mesh->ele_vector.push_back(createElement(elem_type, mat, node_idx_map[index], node_idx_map[index + 1],
															 node_idx_map[index + incHeight]));       // upper left triangle
					mesh->ele_vector.push_back(createElement(elem_type, mat, node_idx_map[index + 1],
															 node_idx_map[index + incHeight + 1],
															 node_idx_map[index + incHeight]));                   // lower right triangle
				}
				if (elem_type == MshElemType::QUAD)
				{
					mesh->ele_vector.push_back(createElement(elem_type, mat, node_idx_map[index], node_idx_map[index + 1],
															 node_idx_map[index + incHeight + 1],
															 node_idx_map[index + incHeight]));
				}
			}
		}

	mesh->ConstructGrid();
	return mesh;
}

MeshLib::CElem* VtkMeshConverter::createElement(MshElemType::type t, int mat, size_t node1, size_t node2, size_t node3, size_t node4)
{
	MeshLib::CElem* elem(new MeshLib::CElem);
	elem->setElementProperties(t);
	elem->SetNodeIndex(0, node1);
	elem->SetNodeIndex(1, node2);
	elem->SetNodeIndex(2, node3);
	if (t ==  MshElemType::QUAD)
		elem->SetNodeIndex(3, node4);
	elem->SetPatchIndex(mat);
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
		elem->getNodeIndices().resize(nElemNodes);

		for (size_t j = 0; j < nElemNodes; j++)
			elem->SetNodeIndex(j, cell->GetPointId(j));

		mesh->ele_vector.push_back(elem);
	}
	mesh->ConstructGrid();
	return mesh;
}
