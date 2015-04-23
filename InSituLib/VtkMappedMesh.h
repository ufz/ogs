/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-27
 * \brief  VtkMappedMesh is a adapter for converting OGS cell connectivity to its
 *         VTK counter part. This generates an incomplete vtkUnstructureGrid
 *         without node coordinates. See VtkMeshNodalCoordinatesTemplate.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKMAPPEDMESH_H_
#define VTKMAPPEDMESH_H_

#include <vtkObject.h>
#include <vtkMappedUnstructuredGrid.h>

#include "MeshEnums.h"

class vtkGenericCell;
namespace MeshLib {
	class Element;
	class Node;
}

namespace InSituLib
{

class VtkMappedMeshImpl : public vtkObject
{
public:
	static VtkMappedMeshImpl *New();
	virtual void PrintSelf(std::ostream &os, vtkIndent indent);
	vtkTypeMacro(VtkMappedMeshImpl, vtkObject)

	void SetNodes(std::vector<MeshLib::Node*> const & nodes);
	void SetElements(std::vector<MeshLib::Element*> const & elements);

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
	VtkMappedMeshImpl();
	~VtkMappedMeshImpl();

private:
	VtkMappedMeshImpl(const VtkMappedMeshImpl &);  // Not implemented.
	void operator=(const VtkMappedMeshImpl &); // Not implemented.

	// const MeshLib::Mesh* _mesh;
	const std::vector<MeshLib::Node*>* _nodes;
	const std::vector<MeshLib::Element*>* _elements;
	vtkIdType NumberOfCells;

	static CellType VtkCellTypeToOGS(int type)
	{
		CellType ogs;
		switch (type)
		{
			case VTK_LINE:
				ogs = CellType::LINE2;
				break;
			case VTK_QUADRATIC_EDGE:
				ogs = CellType::LINE3;
				break;
			case VTK_TRIANGLE:
				ogs = CellType::TRI3;
				break;
			case VTK_QUADRATIC_TRIANGLE:
				ogs = CellType::TRI6;
				break;
			case VTK_QUAD:
				ogs = CellType::QUAD4;
				break;
			case VTK_QUADRATIC_QUAD:
				ogs = CellType::QUAD8;
				break;
			case VTK_BIQUADRATIC_QUAD:
				ogs = CellType::QUAD9;
				break;
			case VTK_HEXAHEDRON:
				ogs = CellType::HEX8;
				break;
			case VTK_QUADRATIC_HEXAHEDRON:
				ogs = CellType::HEX20;
				break;
			case VTK_TRIQUADRATIC_HEXAHEDRON:
				ogs = CellType::HEX27;
				break;
			case VTK_TETRA:
				ogs = CellType::TET4;
				break;
			case VTK_QUADRATIC_TETRA:
				ogs = CellType::TET10;
				break;
			case VTK_WEDGE:
				ogs = CellType::PRISM6;
				break;
			case VTK_QUADRATIC_WEDGE:
				ogs = CellType::PRISM15;
				break;
			case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
				ogs = CellType::PRISM18;
				break;
			case VTK_PYRAMID:
				ogs = CellType::PYRAMID5;
				break;
			case VTK_QUADRATIC_PYRAMID:
				ogs = CellType::PYRAMID13;
				break;
			default:
				ogs = CellType::INVALID;
				break;
		}
		return ogs;
	}
};

vtkMakeMappedUnstructuredGrid(VtkMappedMesh, VtkMappedMeshImpl)

} // end namespace

#endif // VTKMAPPEDMESH_H_
