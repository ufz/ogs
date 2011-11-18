/**
 * \file VtkMeshConverter.h
 * 23/08/2011 KR Initial implementation
 *
 */

#ifndef VTKMESHCONVERTER_H
#define VTKMESHCONVERTER_H

// ** INCLUDES **
#include "msh_mesh.h"

class vtkImageData; // For conversion from Image to QuadMesh
class vtkUnstructuredGrid; // For conversion vom vtk to ogs mesh

namespace MeshLib
{
class CFEMesh;
class CNode;
}

/**
 * \brief Adapter class to convert FEM Mesh to a representation more suited for visualisation purposes
 */
class VtkMeshConverter
{
public:
	/**
	 * Converts greyscale image to a mesh
	 * \param type defines if elements of the new mesh should be triangles or quads (or hexes for 3D)
	 * \param setAsElevation defines how image intensities are interpreted (as elevation or as material IDs)
	 */
	static MeshLib::CFEMesh* convertImgToMesh(vtkImageData* img,
	                                          const std::pair<double,double> &origin,
	                                          const double &scalingFactor,
											  MshElemType::type t,
											  bool setAsElevation);

	/// Converts a vtkUnstructuredGrid object to a CFEMesh
	static MeshLib::CFEMesh* convertUnstructuredGrid(vtkUnstructuredGrid* grid);

private:
	static MeshLib::CElem* createElement(MshElemType::type t, int mat,
		                                 size_t node1, size_t node2, 
										 size_t node3, size_t node4 = 0);
};

#endif // VTKMESHCONVERTER_H
