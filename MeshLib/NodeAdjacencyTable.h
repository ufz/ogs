/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <vector>

#include "Node.h"


namespace MeshLib
{

/// Representation of topological node adjacency.
///
/// The topological sparsity pattern in the context of FEM is defined in terms of
/// supports of the nodal functions. Especially, two nodes i and j are called
/// adjacent if and only if there is a mesh element E including nodes i and j.
/// This information is represented by the NodeAdjacenceTable.
///
/// The topological adjacency of nodes is created by
/// Mesh::setNodesConnectedByElements() which is usually called upon mesh
/// construction.
class
NodeAdjacencyTable
{

public:
    NodeAdjacencyTable() = default;

    explicit
    NodeAdjacencyTable(std::vector<Node*> const& nodes)
    {
        _data.resize(nodes.size());

        createTable(nodes);
    }

    std::size_t size() const
    {
        return _data.size();
    }

    std::size_t getNodeDegree(std::size_t const node_id) const
    {
        return _data[node_id].size();
    }

    std::vector<std::size_t> const& getAdjacentNodes(std::size_t const node_id) const
    {
        return _data[node_id];
    }

    void createTable(std::vector<Node*> const& nodes)
    {
        if (_data.size() != nodes.size())
            _data.resize(nodes.size());

        for (auto n_ptr : nodes)
        {
            std::vector<Node*> const& connected_nodes = n_ptr->getConnectedNodes();
            std::vector<std::size_t>& row = _data[n_ptr->getID()];
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
