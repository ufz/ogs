/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-27
 * \brief  Implementation of the VtkMappedMesh class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMappedMesh.h"

#include "Element.h"

#include <vtkCellType.h>
#include <vtkCellTypes.h>
#include <vtkGenericCell.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkPoints.h>

namespace InSituLib {

vtkStandardNewMacro(VtkMappedMesh)

void VtkMappedMesh::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	// TODO
	os << indent << "NumberOfCells: " << this->NumberOfCells << endl;
}

bool VtkMappedMesh::SetElements(std::vector< MeshLib::Element * > const & elements)
{
	_elements = &elements;

	return true;
}

vtkIdType VtkMappedMesh::GetNumberOfCells()
{
	return this->NumberOfCells;
}

int VtkMappedMesh::GetCellType(vtkIdType cellId)
{
	return (int)(*_elements)[cellId]->getCellType();
}

void VtkMappedMesh::GetCellPoints(vtkIdType cellId, vtkIdList *ptIds)
{
	ptIds->SetNumberOfIds((int)(*_elements)[cellId]->getNNodes());

	std::transform(this->GetElementStart(cellId),
	               this->GetElementEnd(cellId),
	               ptIds->GetPointer(0), NodeToPoint);
}

void VtkMappedMesh::GetPointCells(vtkIdType ptId, *cellIds)
{
	const int targetElement = PointToNode(ptId);
	int *element = this->GetStart();
	int *elementEnd = this->GetEnd();

	cellIds->Reset();

	element = std::find(element, elementEnd, targetElement);
	while (element != elementEnd)
		{
		cellIds->InsertNextId(static_cast<vtkIdType>((element - this->Elements)
		                                             / this->CellSize));
		element = std::find(element, elementEnd, targetElement);
		}
}

int VtkMappedMesh::GetMaxCellSize()
{
	return this->CellSize;
}

void VtkMappedMesh::GetIdsOfCellsOfType(int type,
{
	array->Reset();
	if (type == this->CellType)
		{
		array->SetNumberOfComponents(1);
		array->Allocate(this->NumberOfCells);
		for (vtkIdType i = 0; i < this->NumberOfCells; ++i)
			{
			array->InsertNextValue(i);
			}
		}
}

int VtkMappedMesh::IsHomogeneous()
{
	return 1;
}

void VtkMappedMesh::Allocate(vtkIdType, int)
{
	vtkErrorMacro("Read only container.")
	return;
}

vtkIdType VtkMappedMesh::InsertNextCell(int, vtkIdList*)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

vtkIdType VtkMappedMesh::InsertNextCell(int, vtkIdType, vtkIdType*)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

vtkIdType VtkMappedMesh::InsertNextCell(
		int, vtkIdType, vtkIdType*, vtkIdType, vtkIdType*)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

void VtkMappedMesh::ReplaceCell(vtkIdType, int, vtkIdType*)
{
	vtkErrorMacro("Read only container.")
	return;
}

VtkMappedMesh::VtkMappedMesh()
	: Elements(NULL),
		CellType(VTK_EMPTY_CELL),
		CellSize(0),
		NumberOfCells(0)
{
}

VtkMappedMesh::~VtkMappedMesh()
{
	delete [] this->Elements;
}

} // end namespace
