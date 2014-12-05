/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLIB_NODE_ADJACENCE_TABLE_H_
#define MESHLIB_NODE_ADJACENCE_TABLE_H_

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
    NodeAdjacencyTable(std::vector<Node*> const& nodes, std::size_t hint_adjacency_degree = 0)
    {
        _data.resize(nodes.size());
        if (hint_adjacency_degree > 0)
            for (auto row : _data)
                row.reserve(hint_adjacency_degree);

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

    void createTable(std::vector<Node*> const& nodes, std::size_t hint_adjacency_degree = 0)
    {
        if (_data.size() != nodes.size())
            _data.resize(nodes.size());

        // Allocate temporary space for adjacent nodes.
        std::vector<std::size_t> adjacent_nodes;

        for (auto n_ptr : nodes)
        {
            adjacent_nodes.clear();

            // Get all elements, to which this node is connected.
            std::vector<Element*> const& connected_elements = n_ptr->getElements();

            // And collect all elements' nodes.
            for (auto e : connected_elements)
            {
                Node* const* const single_elem_nodes = e->getNodes();
                std::size_t const nnodes = e->getNNodes();
                for (std::size_t n = 0; n < nnodes; n++)
                    adjacent_nodes.push_back(single_elem_nodes[n]->getID());
            }

            // Copy only unique node ids.
            std::sort(adjacent_nodes.begin(), adjacent_nodes.end());

            std::size_t const node_id = n_ptr->getID();
            std::vector<std::size_t>& row = _data[node_id];
            row.reserve(hint_adjacency_degree);
            std::unique_copy(adjacent_nodes.begin(), adjacent_nodes.end(),
                std::back_inserter(row));
        }
    }

private:
    std::vector<std::vector<std::size_t>> _data;

};

}   // namespace MeshLib

#endif  //MESHLIB_NODE_ADJACENCE_TABLE_H_
