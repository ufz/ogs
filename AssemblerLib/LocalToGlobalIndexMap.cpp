/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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



template <typename ElementIterator>
void
LocalToGlobalIndexMap::findGlobalIndices(
    ElementIterator first, ElementIterator last,
    std::size_t const mesh_id,
    const unsigned comp_id, const unsigned comp_id_write)
{
    _rows.resize(std::distance(first, last), _mesh_subsets.size());

    // For each element find the global indices for node/element
    // components.
    std::size_t elem_id = 0;
    for (ElementIterator e = first; e != last; ++e, ++elem_id)
    {
        std::size_t const nnodes = (*e)->getNNodes();

        LineIndex indices;
        indices.reserve(nnodes);

        for (unsigned n = 0; n < nnodes; n++)
        {
            MeshLib::Location l(mesh_id,
                                MeshLib::MeshItemType::Node,
                                (*e)->getNode(n)->getID());
            indices.push_back(_mesh_component_map.getGlobalIndex(l, comp_id));
        }

        _rows(elem_id, comp_id_write) = std::move(indices);
    }
}


LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
    AssemblerLib::ComponentOrder const order)
    : _mesh_subsets(mesh_subsets), _mesh_component_map(_mesh_subsets, order)
{
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.

    unsigned comp_id = 0;
    for (MeshLib::MeshSubsets const* const mss : _mesh_subsets)
    {
        for (MeshLib::MeshSubset const* const ms : *mss)
        {
            std::size_t const mesh_id = ms->getMeshID();

            findGlobalIndices(ms->elementsBegin(), ms->elementsEnd(), mesh_id,
                              comp_id, comp_id);
        }
        ++comp_id;
    }
}

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubsets*>&& mesh_subsets,
    std::vector<std::size_t> const& original_indices,
    std::vector<MeshLib::Element*> const& elements,
    AssemblerLib::MeshComponentMap&& mesh_component_map)
    : _mesh_subsets(std::move(mesh_subsets)),
      _mesh_component_map(std::move(mesh_component_map))
{
    assert(original_indices.size() == _mesh_subsets.size());
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.

    unsigned comp_id = 0;
    for (MeshLib::MeshSubsets const* const mss : _mesh_subsets)
    {
        if (! mss) continue;
        for (MeshLib::MeshSubset const* const ms : *mss)
        {
            std::size_t const mesh_id = ms->getMeshID();

            findGlobalIndices(elements.cbegin(), elements.cend(), mesh_id,
                              original_indices[comp_id], comp_id);
        }
        ++comp_id;
    }
}

LocalToGlobalIndexMap*
LocalToGlobalIndexMap::deriveBoundaryConstrainedMap(std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
    std::vector<MeshLib::Element*> const& elements) const
{
    DBUG("Construct reduced local to global index map.");

    std::vector<MeshLib::MeshSubsets*> subsets;
    std::vector<std::size_t> orig_idcs;
    unsigned i=0;
    for (auto m : mesh_subsets)
    {
        if (m != nullptr) {
            subsets.push_back(m);
            orig_idcs.push_back(i);
        }
        ++i;
    }

    return new LocalToGlobalIndexMap(std::move(subsets), orig_idcs, elements,
        _mesh_component_map.getSubset(mesh_subsets));
}

std::size_t
LocalToGlobalIndexMap::dofSize() const
{
    return _mesh_component_map.size();
}

std::size_t
LocalToGlobalIndexMap::size() const
{
    return _rows.rows();
}

LocalToGlobalIndexMap::RowColumnIndices
LocalToGlobalIndexMap::operator()(std::size_t const mesh_item_id, const unsigned component_id) const
{
    return RowColumnIndices(_rows(mesh_item_id, component_id),
                            _columns(mesh_item_id, component_id));
}

std::size_t
LocalToGlobalIndexMap::getNumElementDOF(std::size_t const mesh_item_id) const
{
    std::size_t ndof = 0;

    for (unsigned c=0; c<_rows.cols(); ++c)
    {
        ndof += _rows(mesh_item_id, c).size();
    }

    return ndof;
}

#ifndef NDEBUG
std::ostream& operator<<(std::ostream& os, LocalToGlobalIndexMap const& map)
{
    std::size_t const max_lines = 10;
    std::size_t lines_printed = 0;

    os << "Rows of the local to global index map; " << map._rows.size()
        << " rows\n";
    for (std::size_t e=0; e<map.size(); ++e)
    {
        for (std::size_t c=0; c<map.getNumComponents(); ++c)
        {
            auto const& line = map._rows(e, c);

            os << "c" << c << " { ";
            std::copy(line.cbegin(), line.cend(),
                std::ostream_iterator<std::size_t>(os, " "));
            os << " }\n";
        }

        if (lines_printed++ > max_lines)
        {
            os << "...\n";
            break;
        }
    }
    lines_printed = 0;

    os << "Mesh component map:\n" << map._mesh_component_map;
    return os;
}
#endif  // NDEBUG

}   // namespace AssemblerLib
