/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-27
 * \brief  Implementation of the VtkMappedMesh class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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
	switch ((*_elements)[cellId]->getCellType())
	{
		case CellType::INVALID:
			break;
		case CellType::LINE2:
			type = VTK_LINE;
			break;
		case CellType::LINE3:
			type = VTK_QUADRATIC_EDGE;
			break;
		case CellType::TRI3:
			type = VTK_TRIANGLE;
			break;
		case CellType::TRI6:
			type = VTK_QUADRATIC_TRIANGLE;
			break;
		case CellType::QUAD4:
			type = VTK_QUAD;
			break;
		case CellType::QUAD8:
			type = VTK_QUADRATIC_QUAD;
			break;
		case CellType::QUAD9:
			type = VTK_BIQUADRATIC_QUAD;
			break;
		case CellType::HEX8:
			type = VTK_HEXAHEDRON;
			break;
		case CellType::HEX20:
			type = VTK_QUADRATIC_HEXAHEDRON;
			break;
		case CellType::HEX27:
			type = VTK_TRIQUADRATIC_HEXAHEDRON;
			break;
		case CellType::TET4:
			type = VTK_TETRA;
			break;
		case CellType::TET10:
			type = VTK_QUADRATIC_TETRA;
			break;
		case CellType::PRISM6:
			type = VTK_WEDGE;
			break;
		case CellType::PRISM15:
			type = VTK_QUADRATIC_WEDGE;
			break;
		case CellType::PRISM18:
			type = VTK_BIQUADRATIC_QUADRATIC_WEDGE;
			break;
		case CellType::PYRAMID5:
			type = VTK_PYRAMID;
			break;
		case CellType::PYRAMID13:
			type = VTK_QUADRATIC_PYRAMID;
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
		if ((*elem)->getCellType() == VtkCellTypeToOGS(type))
			array->InsertNextValue((*elem)->getID());
	}
}

int VtkMappedMeshImpl::IsHomogeneous()
{
	CellType type = (*(_elements->begin()))->getCellType();
	for (auto elem(_elements->begin()); elem != _elements->end(); ++elem)
		if((*elem)->getCellType() != type)
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
