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
	/// Converts greyscale image to quad mesh
	static MeshLib::CFEMesh* convertImgToMesh(vtkImageData* img, const std::pair<double,double> &origin, const double &scalingFactor);

	/// Converts a vtkUnstructuredGrid object to a CFEMesh
	static MeshLib::CFEMesh* convertUnstructuredGrid(vtkUnstructuredGrid* grid);

private:
	static MeshLib::CElem* createElement(size_t node1, size_t node2, size_t node3);
};

#endif // VTKMESHCONVERTER_H
