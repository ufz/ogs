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
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.
    for (MeshLib::MeshSubsets const* const mss : _mesh_subsets)
    {
        for (MeshLib::MeshSubset const* const ms : *mss)
        {
            std::size_t const mesh_id = ms->getMeshID();

            // For each element find the global indices for node/element
            // components.
            for (auto e = ms->elementsBegin();
                    e != ms->elementsEnd(); ++e)
            {
                std::vector<MeshLib::Location> vec_items;
                std::size_t const nnodes = (*e)->getNNodes();
                vec_items.reserve(nnodes);

                for (std::size_t n = 0; n < nnodes; n++)
                {
                    vec_items.emplace_back(
                        mesh_id,
                        MeshLib::MeshItemType::Node,
                        (*e)->getNode(n)->getID());
                }

                // Save a line of indices for the current element.
                switch (order)
                {
                    case AssemblerLib::ComponentOrder::BY_LOCATION:
                        _rows.push_back(_mesh_component_map.getGlobalIndices<AssemblerLib::ComponentOrder::BY_LOCATION>(vec_items));
                        break;
                    case AssemblerLib::ComponentOrder::BY_COMPONENT:
                        _rows.push_back(_mesh_component_map.getGlobalIndices<AssemblerLib::ComponentOrder::BY_COMPONENT>(vec_items));
                        break;
                }
            }
        }
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
