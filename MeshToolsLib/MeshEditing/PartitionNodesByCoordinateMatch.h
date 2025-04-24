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

#pragma once

#include <optional>
#include <vector>

namespace MeshLib
{
class Node;
}  // namespace MeshLib

namespace GeoLib
{
class AABB;
}  // namespace GeoLib

namespace MeshToolsLib
{
struct NodesPartitionResult
{
    std::vector<MeshLib::Node*> paired_nodes;
    std::optional<std::vector<std::size_t>> id_mapping;
    std::optional<std::vector<MeshLib::Node*>> non_paired_nodes;
};

/**
 * Assuming at least one element in \c node_vector can be paired with an element
 * in \c tool_node_vector based on identical coordinates,
 * this function partitions \c node_vector into two vectors: one containing
 * nodes that can be paired with elements in \c tool_node_vector, and another
 * containing the remaining nodes or empty depending on the argument \c
 * return_non_paired_nodes.
 *
 * @param node_vector      Node vector to be partitioned.
 * @param tool_node_vector Node vector as a tool for partitioning.
 * @param aabb             Axis-aligned bounding box (AABB) for spatial
 *                         filtering and the optional return of non-paired
 *                         nodes.
 * @param return_non_paired_nodes Indicator whether to return unpaired nodes.
 * @return
 *    -# paired nodes of \c node_vector,
 *    -# node ID mapping from the paired node of \c node_vector to its pair in
 *       \c tool_node_vector or nothing if \c return_non_paired_nodes is false,
 *    -# unpaired nodes of \c node_vector or nothing if
 *      \c return_non_paired_nodes is false.
 */
NodesPartitionResult partitionNodesByCoordinateMatch(
    std::vector<MeshLib::Node*> const& node_vector,
    std::vector<MeshLib::Node*> const& tool_node_vector,
    GeoLib::AABB const& aabb,
    bool const return_non_paired_nodes = true);

}  // namespace MeshToolsLib
