/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 24, 2025, 4:39 PM
 */

#include "PartitionNodesByCoordinateMatch.h"

#include "GeoLib/AABB.h"
#include "GeoLib/OctTree.h"
#include "MeshLib/Node.h"

namespace MeshToolsLib
{
NodesPartitionResult partitionNodesByCoordinateMatch(
    std::vector<MeshLib::Node*> const& node_vector,
    std::vector<MeshLib::Node*> const& tool_node_vector,
    GeoLib::AABB const& aabb,
    bool const return_non_paired_nodes)
{
    auto oct_tree = std::unique_ptr<GeoLib::OctTree<MeshLib::Node, 16>>(
        GeoLib::OctTree<MeshLib::Node, 16>::createOctTree(
            aabb.getMinPoint(), aabb.getMaxPoint(), 1e-16));

    // Push all tool nodes into oct_tree:
    for (auto const node : tool_node_vector)
    {
        MeshLib::Node* node_ptr = nullptr;
        oct_tree->addPoint(node, node_ptr);
    }

    // Find the paired nodes in the node_vector
    std::vector<MeshLib::Node*> paired_nodes;
    std::vector<MeshLib::Node*> other_nodes;
    std::vector<std::size_t> ip_mapping;
    for (auto node : node_vector)
    {
        MeshLib::Node* node_ptr = nullptr;
        if (oct_tree->addPoint(node, node_ptr))
        {
            if (return_non_paired_nodes)
            {
                other_nodes.push_back(node);
            }
            continue;
        }
        paired_nodes.push_back(node);
        ip_mapping.push_back(node_ptr->getID());
    }

    if (return_non_paired_nodes)
    {
        return {paired_nodes, ip_mapping, other_nodes};
    }

    return {paired_nodes, std::nullopt, std::nullopt};
}
}  // namespace MeshToolsLib