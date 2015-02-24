/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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

            findGlobalIndices(ms->elementsBegin(), ms->elementsEnd(), mesh_id, order);
        }
    }
}

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
    std::vector<MeshLib::Element*> const& elements,
    AssemblerLib::MeshComponentMap&& mesh_component_map,
    AssemblerLib::ComponentOrder const order)
    : _mesh_subsets(mesh_subsets), _mesh_component_map(std::move(mesh_component_map))
{
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.
    for (MeshLib::MeshSubsets const* const mss : _mesh_subsets)
    {
        for (MeshLib::MeshSubset const* const ms : *mss)
        {
            std::size_t const mesh_id = ms->getMeshID();

            findGlobalIndices(elements.cbegin(), elements.cend(), mesh_id, order);
        }
    }
}

LocalToGlobalIndexMap*
LocalToGlobalIndexMap::deriveBoundaryConstrainedMap(
    std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
    std::vector<MeshLib::Element*> const& elements,
    AssemblerLib::ComponentOrder const order) const
{
    DBUG("Construct reduced local to global index map.");

    return new LocalToGlobalIndexMap(mesh_subsets, elements,
        _mesh_component_map.getSubset(mesh_subsets),
        order);
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

#ifndef NDEBUG
std::ostream& operator<<(std::ostream& os, LocalToGlobalIndexMap const& map)
{
    std::size_t const max_lines = 10;
    std::size_t lines_printed = 0;

    os << "Rows of the local to global index map; " << map._rows.size()
        << " rows\n";
    for (auto line : map._rows)
    {
        if (lines_printed++ > max_lines)
        {
            os << "...\n";
            break;
        }

        os << "{ ";
        std::copy(line.cbegin(), line.cend(),
            std::ostream_iterator<std::size_t>(os, " "));
        os << " }\n";
    }
    lines_printed = 0;

    os << "Mesh component map:\n" << map._mesh_component_map;
    return os;
}
#endif  // NDEBUG

}   // namespace AssemblerLib
