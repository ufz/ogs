/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on October 1, 2021, 2:13 PM
 */

#include "NodePartitionedMesh.h"

namespace MeshLib
{
std::vector<int> getEndNodeIDRanks(
    std::size_t const n_global_nodes,
    std::vector<std::size_t> const& n_regular_base_nodes_at_rank,
    std::vector<std::size_t> const& n_regular_high_order_nodes_at_rank)
{
    std::vector<int> data;

    std::transform(n_regular_base_nodes_at_rank.begin() + 1,
                   n_regular_base_nodes_at_rank.end(),
                   n_regular_high_order_nodes_at_rank.begin() + 1,
                   std::back_inserter(data), std::plus<int>());

    data.push_back(n_global_nodes);

    return data;
}

NodePartitionedMesh::NodePartitionedMesh(
    const std::string& name,
    const std::vector<Node*>& nodes,
    const std::vector<std::size_t>& glb_node_ids,
    const std::vector<Element*>& elements,
    Properties const& properties,
    const std::size_t n_global_base_nodes,
    const std::size_t n_global_nodes,
    const std::size_t n_regular_nodes,
    std::vector<std::size_t>&& n_regular_base_nodes_at_rank,
    std::vector<std::size_t>&& n_regular_high_order_nodes_at_rank)
    : Mesh(name, nodes, elements, true /* compute_element_neighbors */,
           properties),
      _global_node_ids(glb_node_ids),
      _n_global_base_nodes(n_global_base_nodes),
      _n_global_nodes(n_global_nodes),
      _n_regular_nodes(n_regular_nodes),
      _n_regular_base_nodes_at_rank(std::move(n_regular_base_nodes_at_rank)),
      _n_regular_high_order_nodes_at_rank(
          std::move(n_regular_high_order_nodes_at_rank)),
      _end_node_id_at_rank(
          getEndNodeIDRanks(n_global_nodes, _n_regular_base_nodes_at_rank,
                            _n_regular_high_order_nodes_at_rank)),
      _is_single_thread(false)
{
}

bool NodePartitionedMesh::isGhostNode(const std::size_t node_id) const
{
    return node_id >= _n_regular_nodes;
}

std::size_t NodePartitionedMesh::getMaximumNConnectedNodesToNode() const
{
    auto const& nodes_connections =
        MeshLib::calculateNodesConnectedByElements(*this);
    auto const max_connections = std::max_element(
        nodes_connections.cbegin(), nodes_connections.cend(),
        [](auto const& connections_node_a, auto const& connections_node_b)
        { return (connections_node_a.size() < connections_node_b.size()); });
    // Return the number of connected nodes +1 for the node itself.
    return static_cast<std::size_t>(max_connections->size() + 1);
}

std::size_t NodePartitionedMesh::getPartitionID(
    const std::size_t global_node_id) const
{
    return std::upper_bound(std::cbegin(_end_node_id_at_rank),
                            std::cend(_end_node_id_at_rank),
                            global_node_id) -
           _end_node_id_at_rank.begin();
}
}  // namespace MeshLib
