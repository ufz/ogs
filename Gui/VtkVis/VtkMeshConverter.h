/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-23
 * \brief  Definition of the VtkMeshConverter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#ifndef VTKMESHCONVERTER_H
#define VTKMESHCONVERTER_H

namespace MeshLib {
	class Mesh;
}

#include "MeshEnums.h"

class vtkImageData; // For conversion from Image to QuadMesh
class vtkUnstructuredGrid; // For conversion vom vtk to ogs mesh

/// Selection of possible interpretations for intensities
struct UseIntensityAs
{
	enum type {
		ELEVATION,
		MATERIAL,
		NONE
	};
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
										  MshElemType::type elem_type,
										  UseIntensityAs::type intensity_type);

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
										  MshElemType::type elem_type,
										  UseIntensityAs::type intensity_type);

	/// Converts a vtkUnstructuredGrid object to a Mesh
	static MeshLib::Mesh* convertUnstructuredGrid(vtkUnstructuredGrid* grid);

private:
	/// Does the actual mesh generation based on the data given to the public methods.
	static MeshLib::Mesh* constructMesh(const double* pixVal,
									   int* node_idx_map,
									   const bool* visNodes,
									   const double origin[3],
									   const std::size_t &imgHeight,
									   const std::size_t &imgWidth,
									   const double &scalingFactor,
									   MshElemType::type elem_type,
									   UseIntensityAs::type intensity_type);

	static double getExistingValue(const double* img, std::size_t length);
};

#endif // VTKMESHCONVERTER_H
