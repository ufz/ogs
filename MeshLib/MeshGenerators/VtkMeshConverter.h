/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-23
 * \brief  Definition of the VtkMeshConverter class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#ifndef VTKMESHCONVERTER_H
#define VTKMESHCONVERTER_H

#include <boost/optional.hpp>
#include <vtkDataArray.h>
#include <vtkType.h>

#include "logog/include/logog.hpp"

#include "GeoLib/Raster.h"
#include "MeshLib/Node.h"
#include "MeshLib/Location.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

class vtkImageData; // For conversion from Image to QuadMesh
class vtkUnstructuredGrid; // For conversion vom vtk to ogs mesh
class vtkDataArray; // For node/cell properties

namespace MeshLib {

class Mesh;
class Properties;


/**
 * \brief Adapter class to convert FEM Mesh to a representation more suited for visualisation purposes
 */
class VtkMeshConverter
{
public:
	/**
	 * Converts greyscale image to a mesh
	 * \param elem_type defines if elements of the new mesh should be triangles or quads (or hexes for 3D)
	 * \param intensity_type defines how image intensities are interpreted
	 */
	static MeshLib::Mesh* convertImgToMesh(vtkImageData* img,
	                                       const double origin[3],
	                                       const double scalingFactor,
	                                       MeshElemType elem_type,
	                                       UseIntensityAs intensity_type);

	/**
	 * Converts double array with raster values to a mesh
	 * \param elem_type defines if elements of the new mesh should be triangles or quads (or hexes for 3D)
	 * \param intensity_type defines how image intensities are interpreted
	 */
	static MeshLib::Mesh* convertImgToMesh(const double* img,
	                                      GeoLib::RasterHeader const& header,
	                                      MeshElemType elem_type,
	                                      UseIntensityAs intensity_type);

	/// Converts a vtkUnstructuredGrid object to a Mesh
	static MeshLib::Mesh* convertUnstructuredGrid(vtkUnstructuredGrid* grid,
	                                              std::string const& mesh_name = "vtkUnstructuredGrid");

private:
	static MeshLib::Mesh* constructMesh(
		std::vector<double> const& pix_val,
		std::vector<bool> const& pix_vis,
		GeoLib::RasterHeader const& header,
		MeshLib::MeshElemType elem_type,
		MeshLib::UseIntensityAs intensity_type);

	static void convertScalarArrays(vtkUnstructuredGrid &grid, MeshLib::Mesh &mesh);

	static std::vector<MeshLib::Node*> VtkMeshConverter::createNodeVector(
		std::vector<double> const& elevation,
		std::vector<int> &node_idx_map,
		GeoLib::RasterHeader const& header,
		bool use_elevation);

	/// Creates a mesh element vector based on image data
	static std::vector<MeshLib::Element*> createElementVector(
		std::vector<double> const&  pix_val,
		std::vector<bool> const& pix_vis,
		std::vector<MeshLib::Node*> const& nodes,
		std::vector<int> const& node_idx_map,
		std::size_t const imgHeight,
		std::size_t const imgWidth,
		MeshElemType elem_type);

	/// Creates a scalar array/mesh property based on pixel values
	template<typename T>
	static void fillPropertyVector(
		MeshLib::PropertyVector<T> &prop_vec,
		std::vector<double> const& pix_val,
		std::vector<bool> const& pix_vis,
		const std::size_t &imgHeight,
		const std::size_t &imgWidth,
		MeshElemType elem_type)
	{
		for (std::size_t i = 0; i < imgWidth; i++)
			for (std::size_t j = 0; j < imgHeight; j++)
			{
				std::size_t const idx (i*imgHeight+j);
				if (!pix_vis[i*imgHeight+j])
					continue;
				T val (static_cast<T>(pix_val[i*(imgHeight+1)+j]));
				if (elem_type == MeshElemType::TRIANGLE)
				{
					prop_vec.push_back(val);
					prop_vec.push_back(val);
				}
				else if (elem_type == MeshElemType::QUAD)
					prop_vec.push_back(val);
			}
	}

	static void convertArray(vtkDataArray &array,
	                         MeshLib::Properties &properties,
	                         MeshLib::MeshItemType type);

	template<typename T> static void convertTypedArray(vtkDataArray &array,
	                                                   MeshLib::Properties &properties,
	                                                   MeshLib::MeshItemType type)
	{
		vtkIdType const nTuples (array.GetNumberOfTuples());
		int const nComponents (array.GetNumberOfComponents());
		char const*const array_name (array.GetName());

		boost::optional<MeshLib::PropertyVector<T> &> vec
			(properties.createNewPropertyVector<T>(array_name, type, nComponents));
		if (!vec)
		{
			WARN("Array %s could not be converted to PropertyVector.", array_name);
			return;
		}
		vec->reserve(nTuples*nComponents);
		T* data_array = static_cast<T*>(array.GetVoidPointer(0));
		std::copy(&data_array[0], &data_array[nTuples*nComponents], std::back_inserter(*vec));
		return;
	}

	static double getExistingValue(const double* img, std::size_t length);
};

} // end namespace MeshLib

#endif // VTKMESHCONVERTER_H
