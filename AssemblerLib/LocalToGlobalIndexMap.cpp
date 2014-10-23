/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalToGlobalIndexMap.h"

#include "logog/include/logog.hpp"

#include "AssemblerLib/MeshComponentMap.h"
#include "MeshLib/MeshSubsets.h"

namespace AssemblerLib
{

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
    AssemblerLib::ComponentOrder const order)
    : _mesh_subsets(mesh_subsets), _mesh_component_map(_mesh_subsets, order)
{
    switch (order)
    {
        case AssemblerLib::ComponentOrder::BY_LOCATION:
            constructGlobalIndicesForMeshSubsets<AssemblerLib::ComponentOrder::BY_LOCATION>();
        case AssemblerLib::ComponentOrder::BY_COMPONENT:
        default:
            ERR("Unknown type of ComponentOrder passed. Using BY_COMPONENT");
            constructGlobalIndicesForMeshSubsets<AssemblerLib::ComponentOrder::BY_COMPONENT>();
    }
}

std::size_t
LocalToGlobalIndexMap::dofSize() const
{
    return _mesh_component_map.size();
}

std::size_t
LocalToGlobalIndexMap::size() const
{
    return _rows.size();
}

LocalToGlobalIndexMap::RowColumnIndices
LocalToGlobalIndexMap::operator[](std::size_t const mesh_item_id) const
{
    return RowColumnIndices(_rows[mesh_item_id], _columns[mesh_item_id]);
}

LocalToGlobalIndexMap::LineIndex
LocalToGlobalIndexMap::rowIndices(std::size_t const mesh_item_id) const
{
    return _rows[mesh_item_id];
}

LocalToGlobalIndexMap::LineIndex
LocalToGlobalIndexMap::columnIndices(std::size_t const mesh_item_id) const
{
    return _columns[mesh_item_id];
}

}   // namespace AssemblerLib
