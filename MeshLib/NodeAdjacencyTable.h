/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <vector>

#include "Mesh.h"

namespace MeshLib
{
/// Representation of topological node adjacency.
///
/// The topological sparsity pattern in the context of FEM is defined in terms
/// of supports of the nodal functions. Especially, two nodes i and j are called
/// adjacent if and only if there is a mesh element E including nodes i and j.
/// This information is represented by the NodeAdjacenceTable.
///
/// The topological adjacency of nodes is created by
/// MeshLib::calculateNodesConnectedByElements().
class NodeAdjacencyTable final
{
public:
    explicit NodeAdjacencyTable(Mesh const& mesh)
    {
        _data.resize(mesh.getNodes().size());

        createTable(mesh);
    }

    std::size_t size() const
    {
        return _data.size();
    }

    std::size_t getNodeDegree(std::size_t const node_id) const
    {
        return _data[node_id].size();
    }

    std::vector<std::size_t> const& getAdjacentNodes(
        std::size_t const node_id) const
    {
        return _data[node_id];
    }

    void createTable(Mesh const& mesh)
    {
        auto const& nodes = mesh.getNodes();
        if (_data.size() != nodes.size())
        {
            _data.resize(nodes.size());
        }

        auto const& connections =
            MeshLib::calculateNodesConnectedByElements(mesh);
        for (auto const* node : nodes)
        {
            auto const node_id = node->getID();
            auto const& connected_nodes = connections[node_id];
            std::vector<std::size_t>& row = _data[node_id];
            row.reserve(connected_nodes.size());
            std::transform(connected_nodes.cbegin(), connected_nodes.cend(),
                std::back_inserter(row),
                [](Node const* const n) { return n->getID(); });
        }
    }

private:
    std::vector<std::vector<std::size_t>> _data;
};

}   // namespace MeshLib
