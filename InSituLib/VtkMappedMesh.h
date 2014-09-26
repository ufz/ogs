/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-27
 * \brief  Definition of the VtkMappedMesh class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKMAPPEDMESH_H_
#define VTKMAPPEDMESH_H_

#include <vtkObject.h>
#include <vtkMappedUnstructuredGrid.h>

class vtkGenericCell;
namespace MeshLib {
	class Element;
}

namespace InSituLib
{

class VtkMappedMesh : public vtkObject
{
public:
	static VtkMappedMesh *New();
	virtual void PrintSelf(ostream &os, vtkIndent indent);
	vtkTypeMacro(VtkMappedMesh, vtkObject)

	bool SetElements(std::vector< MeshLib::Element * > const & elements);

	// API for vtkMappedUnstructuredGrid's implementation
	vtkIdType GetNumberOfCells();
	int GetCellType(vtkIdType cellId);
	void GetCellPoints(vtkIdType cellId, vtkIdList *ptIds);
	void GetPointCells(vtkIdType ptId, vtkIdList *cellIds);
	int GetMaxCellSize();
	void GetIdsOfCellsOfType(int type, vtkIdTypeArray *array);
	int IsHomogeneous();

	// This container is read only -- these methods do nothing but print a warning.
	void Allocate(vtkIdType numCells, int extSize = 1000);
	vtkIdType InsertNextCell(int type, vtkIdList *ptIds);
	vtkIdType InsertNextCell(int type, vtkIdType npts, vtkIdType *ptIds);
	vtkIdType InsertNextCell(int type, vtkIdType npts, vtkIdType *ptIds,
	                         vtkIdType nfaces, vtkIdType *faces);
	void ReplaceCell(vtkIdType cellId, int npts, vtkIdType *pts);

protected:
	VtkMappedMesh();
	~VtkMappedMesh();

private:
	VtkMappedMesh(const VtkMappedMesh &);  // Not implemented.
	void operator=(const VtkMappedMesh &); // Not implemented.


	const std::vector<MeshLib::Element*>* _elements;
	vtkIdType NumberOfCells;
};

} // end namespace

#endif // VTKMAPPEDMESH_H_
