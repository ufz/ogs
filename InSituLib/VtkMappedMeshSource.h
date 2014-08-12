/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  Definition of the VtkMappedMeshSource class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef _VTKMAPPEDMESHSOURCE
#define _VTKMAPPEDMESHSOURCE

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkNew.h" // For vtkNew
#include <string> // For std::string
#include <vector> // For std::vector

class vtkDataArrayCollection;
class vtkPointData;
class vtkPoints;
namespace MeshLib {
	class Mesh;
}

namespace InSituLib {

class VtkMappedMeshSource : public vtkUnstructuredGridAlgorithm
{
public:
	static VtkMappedMeshSource *New();
	vtkTypeMacro(VtkMappedMeshSource, vtkUnstructuredGridAlgorithm)
	virtual void PrintSelf(ostream &os, vtkIndent indent);

	void SetMesh(const MeshLib::Mesh* mesh) { this->_mesh = mesh; this->Modified(); }
	const MeshLib::Mesh* GetMesh() const { return _mesh; }

protected:
	VtkMappedMeshSource();
	~VtkMappedMeshSource();

	int ProcessRequest(vtkInformation *request, vtkInformationVector **inputVector,
	                   vtkInformationVector *outputVector);
	int RequestData(vtkInformation *, vtkInformationVector **,
	                vtkInformationVector *);
	int RequestInformation(vtkInformation *, vtkInformationVector **,
	                       vtkInformationVector *);

private:
	VtkMappedMeshSource(const VtkMappedMeshSource &); // Not implemented.
	void operator=(const VtkMappedMeshSource &);   // Not implemented.

	const MeshLib::Mesh* _mesh;

	int NumberOfDimensions;
	int NumberOfNodes;
	//int NumberOfElementBlocks;
	std::vector<std::string> NodalVariableNames;
	std::vector<std::string> ElementVariableNames;
	//std::vector<int> ElementBlockIds;

	bool GetCoords();
	vtkNew<vtkPoints> Points;

	bool GetNodalVars();
	vtkNew<vtkPointData> PointData;

	bool GetElems();
	vtkNew<vtkUnstructuredGrid> ElementData;
};

} // Namespace InSituLib

#endif //_VTKMAPPEDMESHSOURCE
