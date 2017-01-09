/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-27
 * \brief  Implementation of the VtkMappedMesh class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMappedMesh.h"

#include <algorithm>

#include <logog/include/logog.hpp>

#include <vtkCellType.h>
#include <vtkCellTypes.h>
#include <vtkGenericCell.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkPoints.h>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"
#include "MeshLib/VtkOGSEnum.h"

namespace MeshLib {

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
    return OGSToVtkCellType((*_elements)[cellId]->getCellType());
}

void VtkMappedMeshImpl::GetCellPoints(vtkIdType cellId, vtkIdList *ptIds)
{
    const MeshLib::Element* const elem = (*_elements)[cellId];
    const unsigned numNodes(elem->getNumberOfNodes());
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
    else if(GetCellType(cellId) == VTK_QUADRATIC_WEDGE)
    {
        std::array<vtkIdType, 15> ogs_nodeIds;
        for (unsigned i=0; i<15; ++i)
            ogs_nodeIds[i] = ptIds->GetId(i);
        for (unsigned i=0; i<3; ++i)
        {
            ptIds->SetId(i, ogs_nodeIds[i+3]);
            ptIds->SetId(i+3, ogs_nodeIds[i]);
        }
        for (unsigned i=0; i<3; ++i)
            ptIds->SetId(6+i, ogs_nodeIds[8-i]);
        for (unsigned i=0; i<3; ++i)
            ptIds->SetId(9+i, ogs_nodeIds[14-i]);
        ptIds->SetId(12, ogs_nodeIds[9]);
        ptIds->SetId(13, ogs_nodeIds[11]);
        ptIds->SetId(14, ogs_nodeIds[10]);
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
        size = std::max(size, (*elem)->getNumberOfNodes());
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
    MeshLib::CellType type = (*(_elements->begin()))->getCellType();
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
