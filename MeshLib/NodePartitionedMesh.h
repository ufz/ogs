/*!
  \file
  \author Wenqing Wang
  \date   2014.06
  \brief  Definition of mesh class for partitioned mesh (by node) for parallel
  computing within the framework of domain decomposition (DDC).

  \copyright
  Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "Mesh.h"
#include "Node.h"

namespace MeshLib
{
/// A subdomain mesh.
class NodePartitionedMesh : public Mesh
{
public:
    // Copy a global mesh for the case of the thread number is one,
    // i.e the global mesh is not partitioned.
    // \param mesh The global mesh
    explicit NodePartitionedMesh(const Mesh& mesh)
        : Mesh(mesh),
          _global_node_ids(mesh.getNumberOfNodes()),
          _n_global_base_nodes(mesh.computeNumberOfBaseNodes()),
          _n_global_nodes(mesh.getNumberOfNodes()),
          _n_regular_nodes(mesh.getNumberOfNodes()),
          _is_single_thread(true)
    {
        for (std::size_t i = 0; i < _nodes.size(); i++)
        {
            _global_node_ids[i] = _nodes[i]->getID();
        }
    }

    /*!
        \brief Constructor
        \param name          Name assigned to the mesh.
        \param nodes         Vector for nodes, which storage looks like:
                             |regular base nodes|ghost base nodes| ...
                              ... regular higher order nodes| ...
                              ... ghost higher order nodes|
        \param glb_node_ids  Global IDs of nodes of a partition.
        \param elements      Vector for elements. Ghost elements are stored
                             after regular (non-ghost) elements.
        \param properties    Mesh property.
        \param n_global_base_nodes Number of the base nodes of the global mesh.
        \param n_global_nodes      Number of all nodes of the global mesh.
        \param n_regular_nodes      Number of all regular nodes.
        \param n_regular_base_nodes_at_rank Numbers of the regular base nodes of
                                           all previous ranks.
        \param n_regular_high_order_nodes_at_rank Numbers of the regular high
                                                 order nodes of all previous
                                                 ranks.
    */
    NodePartitionedMesh(
        const std::string& name,
        const std::vector<Node*>& nodes,
        const std::vector<std::size_t>& glb_node_ids,
        const std::vector<Element*>& elements,
        Properties const& properties,
        const std::size_t n_global_base_nodes,
        const std::size_t n_global_nodes,
        const std::size_t n_regular_nodes,
        std::vector<std::size_t>&& n_regular_base_nodes_at_rank,
        std::vector<std::size_t>&& n_regular_high_order_nodes_at_rank);

    /// Get the number of nodes of the global mesh for linear elements.
    std::size_t getNumberOfGlobalBaseNodes() const
    {
        return _n_global_base_nodes;
    }

    /// Get the number of all nodes of the global mesh.
    std::size_t getNumberOfGlobalNodes() const { return _n_global_nodes; }
    /// Get the global node ID of a node with its local ID.
    std::size_t getGlobalNodeID(const std::size_t node_id) const
    {
        return _global_node_ids[node_id];
    }

    /// Get the number of all regular nodes of the partition.
    std::size_t getNumberOfRegularNodes() const { return _n_regular_nodes; }
    /// Check whether a node with ID of node_id is a ghost node
    bool isGhostNode(const std::size_t node_id) const;

    std::size_t getNumberOfRegularBaseNodesAtRank(int const partition_id) const
    {
        return _n_regular_base_nodes_at_rank[partition_id];
    }

    std::size_t getNumberOfRegularHighOrderNodesAtRank(
        int const partition_id) const
    {
        return _n_regular_high_order_nodes_at_rank[partition_id];
    }

    // TODO I guess that is a simplified version of computeSparsityPattern()
    /// Get the maximum number of connected nodes to node.
    std::size_t getMaximumNConnectedNodesToNode() const;

    std::size_t getPartitionID(const std::size_t global_node_id) const;

    int getNumberOfPartitions() const
    {
        return _n_regular_base_nodes_at_rank.size();
    }

    bool isForSingleThread() const { return _is_single_thread; }

private:
    /// Global IDs of nodes of a partition
    std::vector<std::size_t> _global_node_ids;

    /// Number of the nodes of the global mesh linear interpolations.
    std::size_t _n_global_base_nodes;

    /// Number of all nodes of the global mesh.
    std::size_t _n_global_nodes;

    /// Number of the all regular nodes.
    std::size_t _n_regular_nodes;

    /// Gathered numbers of the regular nodes for linear interpolations of all
    /// partitions.
    std::vector<std::size_t> _n_regular_base_nodes_at_rank;

    /// Gathered numbers of the all regular high order nodes of all partitions.
    std::vector<std::size_t> _n_regular_high_order_nodes_at_rank;

    /// Gathered the end node id of each rank.
    std::vector<int> _end_node_id_at_rank;

    const bool _is_single_thread;
};

}  // namespace MeshLib
