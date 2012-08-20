/**
 * \file VtkMeshConverter.cpp
 * 23/08/2011 KR Initial implementation
 *
 */

#include "VtkMeshConverter.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Edge.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"

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


MeshLib::Mesh* VtkMeshConverter::convertImgToMesh(vtkImageData* img,
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

	MeshLib::Mesh* mesh = constructMesh(pixVal, node_idx_map, visNodes, origin, imgHeight, imgWidth, scalingFactor, elem_type, intensity_type);

	delete [] pixVal;
	delete [] visNodes;
	delete [] node_idx_map;

	return mesh;
}

MeshLib::Mesh* VtkMeshConverter::convertImgToMesh(const double* img,
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

	MeshLib::Mesh* mesh = constructMesh(pixVal, node_idx_map, visNodes, origin, imgHeight, imgWidth, scalingFactor, elem_type, intensity_type);

	delete [] pixVal;
	delete [] visNodes;
	delete [] node_idx_map;

	return mesh;
}

MeshLib::Mesh* VtkMeshConverter::constructMesh(const double* pixVal,
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
	size_t node_idx_count(0);
	const double x_offset(origin[0] - scalingFactor/2.0);
	const double y_offset(origin[1] - scalingFactor/2.0);

	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

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
				MeshLib::Node* node (new MeshLib::Node(x_offset + (scalingFactor * j), y_offset + (scalingFactor * i), zValue));
				nodes.push_back(node);
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
					MeshLib::Tri* tri1 (new MeshLib::Tri(nodes[node_idx_map[index]], nodes[node_idx_map[index + 1]], 
						                                 nodes[node_idx_map[index + incHeight]], mat));	// upper left triangle
					MeshLib::Tri* tri2 (new MeshLib::Tri(nodes[node_idx_map[index + 1]], nodes[node_idx_map[index + incHeight + 1]], 
						                                 nodes[node_idx_map[index + incHeight]], mat));	// lower right triangle
					elements.push_back(tri1);
					elements.push_back(tri2);                 
				}
				if (elem_type == MshElemType::QUAD)
				{
					MeshLib::Quad* quad (new MeshLib::Quad(nodes[node_idx_map[index]], nodes[node_idx_map[index + 1]],
													 nodes[node_idx_map[index + incHeight + 1]], nodes[node_idx_map[index + incHeight]], mat));
					elements.push_back(quad);
				}
			}
		}

	return new MeshLib::Mesh("vtkImageData-Mesh", nodes, elements);
}

MeshLib::Mesh* VtkMeshConverter::convertUnstructuredGrid(vtkUnstructuredGrid* grid)
{
	if (!grid)
		return NULL;

	// set mesh nodes
	const size_t nNodes = grid->GetPoints()->GetNumberOfPoints();
	std::vector<MeshLib::Node*> nodes(nNodes);
	double* coords = NULL;
	for (size_t i = 0; i < nNodes; i++)
	{
		coords = grid->GetPoints()->GetPoint(i);
		nodes[i] = new MeshLib::Node(coords[0], coords[1], coords[2]);
	}

	// set mesh elements
	const size_t nElems = grid->GetNumberOfCells();
	std::vector<MeshLib::Element*> elements(nElems);
	vtkDataArray* scalars = grid->GetCellData()->GetScalars("MaterialIDs");
	for (size_t i = 0; i < nElems; i++)
	{
		MeshLib::Element* elem;
		const size_t nElemNodes (grid->GetCell(i)->GetNumberOfPoints());
		std::vector<unsigned> node_ids(nElemNodes);
		for (size_t j=0; j<nElemNodes; j++)
			node_ids[j] = grid->GetCell(i)->GetPointId(j);
		const unsigned material = (scalars) ? static_cast<int>(scalars->GetComponent(i,0)) : 0;

		int cell_type = grid->GetCellType(i);
		switch (cell_type)
		{
		case VTK_TRIANGLE:      
			elem = new MeshLib::Tri(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], material);
			break;
		case VTK_QUAD:          
			elem = new MeshLib::Quad(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], nodes[node_ids[3]], material);
			break;
		case VTK_TETRA:
			elem = new MeshLib::Quad(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], nodes[node_ids[3]], material);
			break;
		case VTK_HEXAHEDRON:
			elem = new MeshLib::Hex(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], nodes[node_ids[3]], 
				                          nodes[node_ids[4]], nodes[node_ids[5]], nodes[node_ids[6]], nodes[node_ids[7]], material);
			break;
		case VTK_PYRAMID:
			elem = new MeshLib::Pyramid(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]],
				                        nodes[node_ids[3]], nodes[node_ids[4]], material);
			break;
		case VTK_WEDGE:         
			elem = new MeshLib::Prism(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], 
				                      nodes[node_ids[3]], nodes[node_ids[4]], nodes[node_ids[5]], material);
			break;
		default:
			std::cout << "Error in GridAdapter::convertUnstructuredGrid() - Unknown mesh element type \"" << cell_type << "\" ..." << std::endl;
			return NULL;
		}

		elements[i] = elem;
	}

	return new MeshLib::Mesh("vtkUnstructuredGrid", nodes, elements);

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
