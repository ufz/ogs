/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-23
 * \brief  Definition of the VtkMeshConverter class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#ifndef VTKMESHCONVERTER_H
#define VTKMESHCONVERTER_H

#include "MeshLib/MeshEnums.h"
#include "MeshLib/Location.h"

class vtkImageData;
class vtkUnstructuredGrid;
class vtkDataArray;

namespace MeshLib {

class Mesh;
class Properties;

/// Selection of possible interpretations for intensities
enum class UseIntensityAs
{
	ELEVATION,
	MATERIAL,
	NONE
};


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
	                                      const double origin[3],
	                                      const std::size_t imgHeight,
	                                      const std::size_t imgWidth,
	                                      const double &scalingFactor,
	                                      MeshElemType elem_type,
	                                      UseIntensityAs intensity_type);

	/// Converts a vtkUnstructuredGrid object to a Mesh
	static MeshLib::Mesh* convertUnstructuredGrid(vtkUnstructuredGrid* grid, 
	                                              std::string const& mesh_name = "vtkUnstructuredGrid");

private:
	/// Does the actual mesh generation based on the data given to the public methods.
	static MeshLib::Mesh* constructMesh(const double* pixVal,
	                                    int* node_idx_map,
	                                    const bool* visNodes,
	                                    const double origin[3],
	                                    const std::size_t &imgHeight,
	                                    const std::size_t &imgWidth,
	                                    const double &scalingFactor,
	                                    MeshElemType elem_type,
	                                    UseIntensityAs intensity_type);

	/// return the first value that is not a "no data" value (i.e. a value that actually exists in the data)
	static double getExistingValue(const double* img, std::size_t length);

	/// Looks up all existing scalar array in the grid object and converts them (if possible) via convertArray.
	static void convertScalarArrays(vtkUnstructuredGrid &grid, MeshLib::Mesh &mesh);

	/// Converts a vtkDataArray into a std::vector and adds it to mesh properties.
	/// Currently only data arrays containing double or int values (or tuples) are supported.
	static void convertArray(vtkDataArray* array, MeshLib::Properties &properties, MeshLib::MeshItemType type);
};

} // end namespace MeshLib

#endif // VTKMESHCONVERTER_H
