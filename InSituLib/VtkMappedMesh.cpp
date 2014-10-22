/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-27
 * \brief  Implementation of the VtkMappedMesh class.
 *
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMappedMesh.h"

#include <algorithm>

#include "logog/include/logog.hpp"

#include <vtkCellType.h>
#include <vtkCellTypes.h>
#include <vtkGenericCell.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkPoints.h>

#include "Element.h"
#include "MeshLib/Node.h"
#include "MeshEnums.h"

namespace InSituLib {

vtkStandardNewMacro(VtkMappedMesh)
vtkStandardNewMacro(VtkMappedMeshImpl)

void VtkMappedMeshImpl::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "NumberOfCells: " << this->NumberOfCells << endl;
}

void VtkMappedMeshImpl::SetNodes(std::vector<MeshLib::Node*> const & nodes)
{
	_nodes = &nodes;
}

void VtkMappedMeshImpl::SetElements(std::vector<MeshLib::Element*> const & elements)
{
	_elements = &elements;
	this->NumberOfCells = _elements->size();
	this->Modified();
}

vtkIdType VtkMappedMeshImpl::GetNumberOfCells()
{
	return this->NumberOfCells;
}

int VtkMappedMeshImpl::GetCellType(vtkIdType cellId)
{
	int type = 0;
	switch ((*_elements)[cellId]->getGeomType())
	{
		case MeshElemType::INVALID:
			break;
		case MeshElemType::LINE:
			type = VTK_LINE;
			break;
		case MeshElemType::TRIANGLE:
			type = VTK_TRIANGLE;
			break;
		case MeshElemType::QUAD:
			type = VTK_QUAD;
			break;
		case MeshElemType::HEXAHEDRON:
			type = VTK_HEXAHEDRON;
			break;
		case MeshElemType::TETRAHEDRON:
			type = VTK_TETRA;
			break;
		case MeshElemType::PRISM:
			type = VTK_WEDGE;
			break;
		case MeshElemType::PYRAMID:
			type = VTK_PYRAMID;
			break;
	}
	return type;
}

void VtkMappedMeshImpl::GetCellPoints(vtkIdType cellId, vtkIdList *ptIds)
{
	const MeshLib::Element* const elem = (*_elements)[cellId];
	const unsigned numNodes(elem->getNNodes());
	const MeshLib::Node* const* nodes = (*_elements)[cellId]->getNodes();
	ptIds->SetNumberOfIds(numNodes);

	for (unsigned i(0); i < numNodes; ++i)
		ptIds->SetId(i, nodes[i]->getID());

	if(GetCellType(cellId) == VTK_WEDGE)
	{
		for (unsigned i=0; i<3; ++i)
		{
			const unsigned prism_swap_id = ptIds->GetId(i);
			ptIds->SetId(i, ptIds->GetId(i+3));
			ptIds->SetId(i+3, prism_swap_id);
		}
	}
}

void VtkMappedMeshImpl::GetPointCells(vtkIdType ptId, vtkIdList *cellIds)
{
	cellIds->Reset();

	auto elements((*_nodes)[ptId]->getElements());
	for (auto elem(elements.begin()); elem != elements.end(); ++elem)
		cellIds->InsertNextId((*elem)->getID());
}

int VtkMappedMeshImpl::GetMaxCellSize()
{
	unsigned int size = 0;
	for (auto elem(_elements->begin()); elem != _elements->end(); ++elem)
		size = std::max(size, (*elem)->getNNodes());
	return size;
}

void VtkMappedMeshImpl::GetIdsOfCellsOfType(int type, vtkIdTypeArray *array)
{
	array->Reset();

	for (auto elem(_elements->begin()); elem != _elements->end(); ++elem)
	{
		if ((*elem)->getGeomType() == VtkCellTypeToOGS(type))
			array->InsertNextValue((*elem)->getID());
	}
}

int VtkMappedMeshImpl::IsHomogeneous()
{
	MeshElemType type = (*(_elements->begin()))->getGeomType();
	for (auto elem(_elements->begin()); elem != _elements->end(); ++elem)
		if((*elem)->getGeomType() != type)
			return 0;
	return 1;
}

void VtkMappedMeshImpl::Allocate(vtkIdType, int)
{
	vtkErrorMacro("Read only container.")
	return;
}

vtkIdType VtkMappedMeshImpl::InsertNextCell(int, vtkIdList*)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

vtkIdType VtkMappedMeshImpl::InsertNextCell(int, vtkIdType, vtkIdType*)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

vtkIdType VtkMappedMeshImpl::InsertNextCell(
		int, vtkIdType, vtkIdType*, vtkIdType, vtkIdType*)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

void VtkMappedMeshImpl::ReplaceCell(vtkIdType, int, vtkIdType*)
{
	vtkErrorMacro("Read only container.")
	return;
}

VtkMappedMeshImpl::VtkMappedMeshImpl()
	: _nodes(NULL), _elements(NULL), NumberOfCells(0)
{
}

VtkMappedMeshImpl::~VtkMappedMeshImpl()
{
	// delete [] this->Elements;
}

} // end namespace
