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
#include <vtkFloatArray.h>


GridAdapter* VtkMeshConverter::convertImgToMesh(vtkImageData* img,
                                                     const double origin[3],
                                                     const double scalingFactor,
													 MshElemType::type elem_type,
													 UseIntensityAs::type intensity_type)
{
	if ((elem_type != MshElemType::TRIANGLE) && (elem_type != MshElemType::QUAD))
	{
		std::cout << "Error in VtkMeshConverter::convertImgToMesh() - Invalid Mesh Element Type..." << std::endl;
		return NULL;
	}

	vtkSmartPointer<vtkFloatArray> pixelData = vtkSmartPointer<vtkFloatArray>(
	        vtkFloatArray::SafeDownCast(img->GetPointData()->GetScalars()));
	int* dims = img->GetDimensions();

	const size_t imgHeight = dims[0];
	const size_t imgWidth  = dims[1];
	const size_t incHeight = imgHeight+1;
	const size_t incWidth  = imgWidth+1;
	double* pixVal (new double[incHeight * incWidth]);
	bool* visNodes(new bool[incWidth * incHeight]);
	int* node_idx_map(new int[incWidth * incHeight]);

	for (size_t j = 0; j < incHeight; j++)
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
			int nTuple = pixelData->GetNumberOfComponents();
			double* colour;
			if (nTuple == 2)	//Grey+Alpha
			{
				colour = pixelData->GetTuple2(img_idx);
				pixVal[index] = colour[0];
			
			}
			else if (nTuple == 4)	//RGBA
			{
				colour = pixelData->GetTuple4(img_idx);
				pixVal[index] = 0.3 * colour[0] + 0.6 * colour[1] + 0.1 * colour[2];
			}
			else
			{
				std::cout << "Unsupported pixel composition!" << std::endl;
				return NULL;
			}
			
			visNodes[index] = (colour[nTuple-1] > 0);
			node_idx_map[index]=-1;
		}
		pixVal[(i+2)*incHeight-1]=0;
		visNodes[(i+2)*incHeight-1]=false;
		node_idx_map[(i+2)*incHeight-1]=-1;
	}
	
	GridAdapter* mesh = constructMesh(pixVal, node_idx_map, visNodes, origin, imgHeight, imgWidth, scalingFactor, elem_type, intensity_type);

	delete [] pixVal;
	delete [] visNodes;
	delete [] node_idx_map;

	return mesh;
}

GridAdapter* VtkMeshConverter::convertImgToMesh(const double* img,
													 const double origin[3],
													 const size_t imgHeight,
													 const size_t imgWidth,
													 const double &scalingFactor,
													 MshElemType::type elem_type,
													UseIntensityAs::type intensity_type)
{
	const size_t incHeight = imgHeight+1;
	const size_t incWidth  = imgWidth+1;
	double* pixVal (new double[incHeight * incWidth]);
	bool* visNodes(new bool[incWidth * incHeight]);
	int* node_idx_map(new int[incWidth * incHeight]);

	double noDataValue = getExistingValue(img, imgWidth*imgHeight);

	for (size_t j = 0; j < imgHeight; j++)
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
			if (img[img_idx] == -9999)
			{
				visNodes[index] = false;
				pixVal[index] = noDataValue;
			}
			else
			{
				pixVal[index] = img[img_idx];
				visNodes[index] = true;
			}

			node_idx_map[index]=-1;
		}
		pixVal[(i+2)*incHeight-1]=0;
		visNodes[(i+2)*incHeight-1]=false;
		node_idx_map[(i+2)*incHeight-1]=-1;
	}
	
	GridAdapter* mesh = constructMesh(pixVal, node_idx_map, visNodes, origin, imgHeight, imgWidth, scalingFactor, elem_type, intensity_type);

	delete [] pixVal;
	delete [] visNodes;
	delete [] node_idx_map;

	return mesh;
}

GridAdapter* VtkMeshConverter::constructMesh(const double* pixVal,
												  int* node_idx_map,
												  const bool* visNodes,
												  const double origin[3],
                                                  const size_t &imgHeight,
												  const size_t &imgWidth,
                                                  const double &scalingFactor,
										 		  MshElemType::type elem_type,
												  UseIntensityAs::type intensity_type)
{
	const size_t incHeight = imgHeight+1;
	const size_t incWidth  = imgWidth+1;
	GridAdapter* grid = new GridAdapter();
	size_t node_idx_count(0);
	const double x_offset(origin[0] - scalingFactor/2.0);
	const double y_offset(origin[1] - scalingFactor/2.0);

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
				grid->addNode(new GEOLIB::Point(x_offset + (scalingFactor * j),
										        y_offset + (scalingFactor * i),
										        zValue));
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
					grid->addElement(createElement(elem_type, mat, node_idx_map[index], node_idx_map[index + 1],
													 node_idx_map[index + incHeight]));       // upper left triangle
					grid->addElement(createElement(elem_type, mat, node_idx_map[index + 1],
													 node_idx_map[index + incHeight + 1],
													 node_idx_map[index + incHeight]));                   // lower right triangle
				}
				if (elem_type == MshElemType::QUAD)
				{
					grid->addElement(createElement(elem_type, mat, node_idx_map[index], node_idx_map[index + 1],
													 node_idx_map[index + incHeight + 1],
													 node_idx_map[index + incHeight]));
				}
			}
		}

	return grid;
}

GridAdapter::Element* VtkMeshConverter::createElement(MshElemType::type t, int mat, size_t node1, size_t node2, size_t node3, size_t node4)
{
	GridAdapter::Element* elem = new GridAdapter::Element();
	elem->material = mat;
	elem->type = t;
	std::vector<size_t> nodes;
	nodes.push_back(node1);
	nodes.push_back(node2);
	nodes.push_back(node3);
	if (t ==  MshElemType::QUAD)
		nodes.push_back(node4);
	elem->nodes = nodes;
	return elem;
}

GridAdapter* VtkMeshConverter::convertUnstructuredGrid(vtkUnstructuredGrid* grid)
{
	if (!grid)
		return NULL;

	GridAdapter* mesh = new GridAdapter();

	const size_t nNodes = grid->GetPoints()->GetNumberOfPoints();
	const size_t nElems = grid->GetNumberOfCells();

	// set mesh nodes
	double* coords = NULL;
	for (size_t i = 0; i < nNodes; i++)
	{
		coords = grid->GetPoints()->GetPoint(i);
		mesh->addNode(new GEOLIB::Point(coords[0], coords[1], coords[2]));
	}

	// set mesh elements
	vtkCell* cell(NULL);
	vtkDataArray* scalars = grid->GetCellData()->GetScalars("MaterialIDs");
	for (size_t i = 0; i < nElems; i++)
	{
		GridAdapter::Element* elem = new GridAdapter::Element();

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
		case VTK_PYRAMID:       elem_type = MshElemType::PYRAMID;
			break;
		}

		if (elem_type != MshElemType::INVALID)
		{
			elem->type = elem_type;
			if (scalars)
				elem->material = static_cast<int>(scalars->GetComponent(i,0));
		}
		else
		{
			std::cout << "Error in GridAdapter::convertUnstructuredGrid() - Unknown mesh element type \"" << cell_type << "\" ..." << std::endl;
			return NULL;
		}

		cell = grid->GetCell(i);
		size_t nElemNodes = cell->GetNumberOfPoints();
		std::vector<size_t> nodes;
		for (size_t j = 0; j < nElemNodes; j++)
			nodes.push_back(cell->GetPointId(j));

		elem->nodes = nodes;
		mesh->addElement(elem);
	}
	return mesh;
}

double VtkMeshConverter::getExistingValue(const double* img, size_t length)
{
	for (size_t i=0; i<length; i++)
	{
		if (img[i] != -9999)
			return img[i];
	}
	return -9999;
}
