/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-23
 * \brief  Implementation of the VtkMeshConverter class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMeshConverter.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Prism.h"

// Conversion from Image to QuadMesh
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

// Conversion from vtkUnstructuredGrid
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>

namespace MeshLib {

MeshLib::Mesh* VtkMeshConverter::convertImgToMesh(vtkImageData* img,
                                                  const double origin[3],
                                                  const double scalingFactor,
                                                  MeshElemType elem_type,
                                                  UseIntensityAs intensity_type)
{
	if ((elem_type != MeshElemType::TRIANGLE) && (elem_type != MeshElemType::QUAD))
	{
		ERR("VtkMeshConverter::convertImgToMesh(): Invalid Mesh Element Type.");
		return nullptr;
	}

	vtkSmartPointer<vtkDataArray> pixelData = vtkSmartPointer<vtkDataArray>(img->GetPointData()->GetScalars());
	int* dims = img->GetDimensions();
	int nTuple = pixelData->GetNumberOfComponents();
	if (nTuple < 1 || nTuple > 4)
	{
		ERR("VtkMeshConverter::convertImgToMesh(): Unsupported pixel composition!");
		return nullptr;
	}

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
			double* colour = pixelData->GetTuple(img_idx);
			if (nTuple < 3)	// Grey (+ Alpha)
				pixVal[index] = colour[0];
			else			// RGB(A)
				pixVal[index] = 0.3 * colour[0] + 0.6 * colour[1] + 0.1 * colour[2];

			// is current pixel visible
			if (nTuple == 2 || nTuple == 4)
				visNodes[index] = (colour[nTuple-1] != 0);
			else
				visNodes[index] = true;

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
                                                  MeshElemType elem_type,
                                                  UseIntensityAs intensity_type)
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
                                               MeshElemType elem_type,
                                               UseIntensityAs intensity_type)
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
				double zValue = (intensity_type == UseIntensityAs::ELEVATION) ? pixVal[index] : origin[2];
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
				if (elem_type == MeshElemType::TRIANGLE)
				{
					MeshLib::Node** tri1_nodes = new MeshLib::Node*[3];
					tri1_nodes[0] = nodes[node_idx_map[index]];
					tri1_nodes[1] = nodes[node_idx_map[index+1]];
					tri1_nodes[2] = nodes[node_idx_map[index+incHeight]];

					MeshLib::Node** tri2_nodes = new MeshLib::Node*[3];
					tri2_nodes[0] = nodes[node_idx_map[index+1]];
					tri2_nodes[1] = nodes[node_idx_map[index+incHeight+1]];
					tri2_nodes[2] = nodes[node_idx_map[index+incHeight]];

					elements.push_back(new MeshLib::Tri(tri1_nodes, mat)); // upper left triangle
					elements.push_back(new MeshLib::Tri(tri2_nodes, mat)); // lower right triangle
				}
				if (elem_type == MeshElemType::QUAD)
				{
					MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
					quad_nodes[0] = nodes[node_idx_map[index]];
					quad_nodes[1] = nodes[node_idx_map[index + 1]];
					quad_nodes[2] = nodes[node_idx_map[index + incHeight + 1]];
					quad_nodes[3] = nodes[node_idx_map[index + incHeight]];
					elements.push_back(new MeshLib::Quad(quad_nodes, mat));
				}
			}
		}

	if (elements.empty())
		return nullptr;

	return new MeshLib::Mesh("RasterDataMesh", nodes, elements); // the name is only a temp-name, the name given in the dialog is set later
}

MeshLib::Mesh* VtkMeshConverter::convertUnstructuredGrid(vtkUnstructuredGrid* grid, std::string const& mesh_name)
{
	if (!grid)
		return nullptr;

	// set mesh nodes
	const size_t nNodes = grid->GetPoints()->GetNumberOfPoints();
	std::vector<MeshLib::Node*> nodes(nNodes);
	double* coords = nullptr;
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
		case VTK_LINE: {
			MeshLib::Node** line_nodes = new MeshLib::Node*[2];
			line_nodes[0] = nodes[node_ids[0]];
			line_nodes[1] = nodes[node_ids[1]];
			elem = new MeshLib::Line(line_nodes, material);
			break;
		}
		case VTK_TRIANGLE: {
			MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
			for (unsigned k(0); k<3; k++)
				tri_nodes[k] = nodes[node_ids[k]];
			elem = new MeshLib::Tri(tri_nodes, material);
			break;
		}
		case VTK_QUAD: {
			MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
			for (unsigned k(0); k<4; k++)
				quad_nodes[k] = nodes[node_ids[k]];
			elem = new MeshLib::Quad(quad_nodes, material);
			break;
		}
		case VTK_PIXEL: {
			MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
			quad_nodes[0] = nodes[node_ids[0]];
			quad_nodes[1] = nodes[node_ids[1]];
			quad_nodes[2] = nodes[node_ids[3]];
			quad_nodes[3] = nodes[node_ids[2]];
			elem = new MeshLib::Quad(quad_nodes, material);
			break;
		}
		case VTK_TETRA: {
			MeshLib::Node** tet_nodes = new MeshLib::Node*[4];
			for (unsigned k(0); k<4; k++)
				tet_nodes[k] = nodes[node_ids[k]];
			elem = new MeshLib::Tet(tet_nodes, material);
			break;
		}
		case VTK_HEXAHEDRON: {
			MeshLib::Node** hex_nodes = new MeshLib::Node*[8];
			for (unsigned k(0); k<8; k++)
				hex_nodes[k] = nodes[node_ids[k]];
			elem = new MeshLib::Hex(hex_nodes, material);
			break;
		}
		case VTK_VOXEL: {
			MeshLib::Node** voxel_nodes = new MeshLib::Node*[8];
			voxel_nodes[0] = nodes[node_ids[0]];
			voxel_nodes[1] = nodes[node_ids[1]];
			voxel_nodes[2] = nodes[node_ids[3]];
			voxel_nodes[3] = nodes[node_ids[2]];
			voxel_nodes[4] = nodes[node_ids[4]];
			voxel_nodes[5] = nodes[node_ids[5]];
			voxel_nodes[6] = nodes[node_ids[7]];
			voxel_nodes[7] = nodes[node_ids[6]];
			elem = new MeshLib::Hex(voxel_nodes, material);
			break;
		}
		case VTK_PYRAMID: {
			MeshLib::Node** pyramid_nodes = new MeshLib::Node*[5];
			for (unsigned k(0); k<5; k++)
				pyramid_nodes[k] = nodes[node_ids[k]];
			elem = new MeshLib::Pyramid(pyramid_nodes, material);
			break;
		}
		case VTK_WEDGE: {
			MeshLib::Node** prism_nodes = new MeshLib::Node*[6];
			for (unsigned k(0); k<6; k++)
				prism_nodes[k] = nodes[node_ids[k]];
			elem = new MeshLib::Prism(prism_nodes, material);
			break;
		}
		default:
			ERR("VtkMeshConverter::convertUnstructuredGrid(): Unknown mesh element type \"%d\".", cell_type);
			return nullptr;
		}

		elements[i] = elem;
	}

	MeshLib::Mesh* mesh = new MeshLib::Mesh(mesh_name, nodes, elements);
	convertScalarArrays(*grid, *mesh);

	return mesh;
}

void VtkMeshConverter::convertScalarArrays(vtkUnstructuredGrid &grid, MeshLib::Mesh &mesh)
{
	vtkPointData* point_data = grid.GetPointData();
	unsigned const n_point_arrays = static_cast<unsigned>(point_data->GetNumberOfArrays());
	for (unsigned i=0; i<n_point_arrays; ++i)
		convertArray(point_data->GetArray(i), mesh.getProperties(), MeshLib::MeshItemType::Node);

	vtkCellData* cell_data = grid.GetCellData();
	unsigned const n_cell_arrays = static_cast<unsigned>(cell_data->GetNumberOfArrays());
	for (unsigned i=0; i<n_cell_arrays; ++i)
		convertArray(cell_data->GetArray(i), mesh.getProperties(), MeshLib::MeshItemType::Cell);
}

void VtkMeshConverter::convertArray(vtkDataArray* array, MeshLib::Properties &properties, MeshLib::MeshItemType type)
{
	vtkIdType const nTuples (array->GetNumberOfTuples());
	int const nComponents (array->GetNumberOfComponents());
	char const*const array_name (array->GetName());

	vtkDoubleArray* double_array = vtkDoubleArray::SafeDownCast(array);
	if (double_array)
	{
		boost::optional<MeshLib::PropertyVector<double> &> vec
			(properties.createNewPropertyVector<double>(array_name, type, nComponents));
		if (!vec)
			return;
		vec->reserve(nTuples*nComponents);
		double* data_array = static_cast<double*>(double_array->GetVoidPointer(0));
		std::copy(&data_array[0], &data_array[nTuples*nComponents], std::back_inserter(*vec));
		return;
	}

	vtkIntArray* int_array = vtkIntArray::SafeDownCast(array);
	if (int_array)
	{
		boost::optional<MeshLib::PropertyVector<int> &> vec
			(properties.createNewPropertyVector<int>(array_name, type, nComponents));
		if (!vec)
			return;
		vec->reserve(nTuples*nComponents);
		int* data_array = static_cast<int*>(int_array->GetVoidPointer(0));
		std::copy(&data_array[0], &data_array[nTuples*nComponents], std::back_inserter(*vec));
		return;
	}

	ERR ("Array \"%s\" in VTU file uses unsupported data type.", array->GetName());
	return;
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

} // end namespace MeshLib
