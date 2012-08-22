/**
 * \file VtkMeshSource.h
 * 19/03/2010 KR Initial implementation
 *
 */

#ifndef VTKMESHSOURCE_H
#define VTKMESHSOURCE_H

// ** INCLUDES **
#include <map>

#include "VtkAlgorithmProperties.h"
#include <vtkUnstructuredGridAlgorithm.h>

class VtkColorLookupTable;

namespace MeshLib {
	class Mesh;
}

/**
 * \brief VTK source object for the visualisation of unstructured grids
 */
class VtkMeshSource : public vtkUnstructuredGridAlgorithm, public VtkAlgorithmProperties
{
public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkMeshSource* New();

	vtkTypeRevisionMacro(VtkMeshSource, vtkUnstructuredGridAlgorithm);

	const char* GetMaterialArrayName() const { return _matName; }

	/// Returns the base object of this grid
	const MeshLib::Mesh* GetMesh() { return this->_grid; }

	/// Sets the grid object that should be visualized
	void SetMesh(const MeshLib::Mesh* grid) { _grid = grid; }

	/// Prints the mesh data to an output stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	/**
	 * \brief Generates random colors based on the material scalar value.
	 * Each element of the mesh is assigned an RGB-value based on its material group.
	 * This method should only be called after setMesh()!
	 */
	//ogsUserPropertyMacro(ColorByMaterial,bool);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkMeshSource();
	~VtkMeshSource();

	/// Computes the unstructured grid data object.
	int RequestData(vtkInformation* request,
	                vtkInformationVector** inputVector,
	                vtkInformationVector* outputVector);

	const MeshLib::Mesh* _grid;

private:
	const char* _matName;
};

#endif // VTKMESHSOURCE_H
