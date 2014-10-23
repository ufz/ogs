/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
#define ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_

#include <vector>

#include "logog/include/logog.hpp"

#include "AssemblerLib/MeshComponentMap.h"
#include "MathLib/LinAlg/RowColumnIndices.h"
#include "MeshLib/MeshSubsets.h"

namespace AssemblerLib
{

/// Row and column indices in global linear algebra objects for each mesh item.
///
/// The row and column indices in global matrix and rhs vector for example,
/// required for addition of local contributions from every mesh item (like node
/// or cell) to global objects.
///
/// The number of rows should be equal to the number of mesh items and the
/// number of columns should be equal to the number of the components on that
/// mesh item.
class LocalToGlobalIndexMap
{
public:
    typedef MathLib::RowColumnIndices<std::size_t> RowColumnIndices;
    typedef RowColumnIndices::LineIndex LineIndex;

public:
    /* \todo Extend the constructor for parallel meshes.
    LocalToGlobalIndexMap(
        std::vector<LineIndex> const& rows,
        std::vector<LineIndex> const & columns)
        : _rows(rows), _columns(columns)
    {
        assert(_rows.size() == _columns.size());
        assert(!_rows.empty());
    }
    */

    /// Creates a MeshComponentMap internally and stores the global indices for
    /// each mesh element of the given mesh_subsets.
    explicit LocalToGlobalIndexMap(
        std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
        AssemblerLib::ComponentOrder const order =
            AssemblerLib::ComponentOrder::BY_COMPONENT)
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
                    switch (order)
                    {
                        case AssemblerLib::ComponentOrder::BY_COMPONENT:
                            appendGlobalIndices<
                                AssemblerLib::ComponentOrder::BY_COMPONENT>(
                                    mesh_id, **e);
                              break;
                        case AssemblerLib::ComponentOrder::BY_LOCATION:
                            appendGlobalIndices<
                                AssemblerLib::ComponentOrder::BY_LOCATION>(
                                    mesh_id, **e);
                              break;
                        default:
                            ERR("Unknown type of ComponentOrder passed.");
                    }
                }
            }
        }
    }

    /// Returns total number of degrees of freedom.
    std::size_t dofSize() const
    {
        return _mesh_component_map.size();
    }

    std::size_t size() const
    {
        return _rows.size();
    }

    RowColumnIndices operator[](std::size_t const mesh_item_id) const
    {
        return RowColumnIndices(_rows[mesh_item_id], _columns[mesh_item_id]);
    }

    LineIndex rowIndices(std::size_t const mesh_item_id) const
    {
        return _rows[mesh_item_id];
    }

    LineIndex columnIndices(std::size_t const mesh_item_id) const
    {
        return _columns[mesh_item_id];
    }

private:
    /// Append global indices for the element e in given component order Order.
    template <AssemblerLib::ComponentOrder Order>
    void appendGlobalIndices(std::size_t const mesh_id, MeshLib::Element const& e)
    {
        std::vector<MeshLib::Location> vec_items;
        std::size_t const nnodes = e.getNNodes();
        vec_items.reserve(nnodes);

        for (std::size_t n = 0; n < nnodes; n++)
        {
            vec_items.emplace_back(
                mesh_id,
                MeshLib::MeshItemType::Node,
                e.getNode(n)->getID());
        }

        // Save a line of indices for the current element.
        _rows.push_back(_mesh_component_map.getGlobalIndices<Order>(vec_items));
    }

private:
    std::vector<MeshLib::MeshSubsets*> const& _mesh_subsets;
    AssemblerLib::MeshComponentMap _mesh_component_map;

    // _rows contains for each element vector of global indices to
    // node/element process variables.
    std::vector<LineIndex> _rows;
    std::vector<LineIndex> const& _columns = _rows;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
