/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  VtkMappedMeshSource is a souce class to transform OGS meshes into complete
 *         vtkUnstructuredGrids.
 * Usage:
 * \code
 * vtkNew<InSituLib::VtkMappedMeshSource> vtkSource;
 * vtkSource->SetMesh(mesh);
 * vtkSource->Update();
 * vtkUnstructuredGrid* output = vtkSource->GetOutput();
 * \endcode
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef _VTKMAPPEDMESHSOURCE
#define _VTKMAPPEDMESHSOURCE

#include <string>
#include <vector>

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkNew.h"

class vtkCellData;
class vtkDataArrayCollection;
class vtkPointData;
class vtkPoints;
namespace MeshLib {
	class Mesh;
	class Properties;
}

namespace InSituLib {

class VtkMappedMeshSource : public vtkUnstructuredGridAlgorithm
{
public:
	static VtkMappedMeshSource *New();
	vtkTypeMacro(VtkMappedMeshSource, vtkUnstructuredGridAlgorithm)
	virtual void PrintSelf(std::ostream &os, vtkIndent indent);

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
	void operator=(const VtkMappedMeshSource &);      // Not implemented.

	template<typename T> bool addProperty(MeshLib::Properties const& properties,
	                                      std::string const& prop_name) const;

	const MeshLib::Mesh* _mesh;

	int NumberOfDimensions;
	int NumberOfNodes;
	std::vector<std::string> NodalVariableNames;
	std::vector<std::string> ElementVariableNames;

	vtkNew<vtkPoints> Points;
	vtkNew<vtkPointData> PointData;
	vtkNew<vtkCellData> CellData;
};

} // Namespace InSituLib

#endif //_VTKMAPPEDMESHSOURCE
