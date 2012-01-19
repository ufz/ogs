/**
 * \file VtkMeshConverter.h
 * 23/08/2011 KR Initial implementation
 *
 */

#ifndef VTKMESHCONVERTER_H
#define VTKMESHCONVERTER_H

//#include <utility>
//#include "MSHEnums.h"
#include "GridAdapter.h"

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
	 * \parelem_type defines if elements of the new mesh should be triangles or quads (or hexes for 3D)
	 * \param intensity_type defines how image intensities are interpreted
	 */
	static GridAdapter* convertImgToMesh(vtkImageData* img,
									      const std::pair<double,double> &origin,
	                                      const double &scalingFactor,
										  MshElemType::type elem_type,
										  UseIntensityAs::type intensity_type);

	/**
	 * Converts double array with raster values to a mesh
	 * \parelem_type defines if elements of the new mesh should be triangles or quads (or hexes for 3D)
	 * \param intensity_type defines how image intensities are interpreted
	 */
	static GridAdapter* convertImgToMesh(const double* img,
	                                      const std::pair<double,double> &origin,
										  const size_t imgHeight,
										  const size_t imgWidth,
	                                      const double &scalingFactor,
										  MshElemType::type elem_type,
										  UseIntensityAs::type intensity_type);

	/// Converts a vtkUnstructuredGrid object to a CFEMesh
	static GridAdapter* convertUnstructuredGrid(vtkUnstructuredGrid* grid);

private:
	/// Does the actual mesh generation based on the data given to the public methods.
	static GridAdapter* constructMesh(const double* pixVal,
									   int* node_idx_map,
									   const bool* visNodes,
									   const std::pair<double,double> &origin,
									   const size_t &imgHeight,
									   const size_t &imgWidth,
									   const double &scalingFactor,
									   MshElemType::type elem_type,
									   UseIntensityAs::type intensity_type);

	/// Creates a mesh element based on the given data.
	static GridAdapter::Element* createElement(MshElemType::type t, int mat,
		                                 size_t node1, size_t node2, 
										 size_t node3, size_t node4 = 0);

	static double getExistingValue(const double* img, size_t length);
};

#endif // VTKMESHCONVERTER_H
