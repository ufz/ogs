/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>

#pragma once

namespace MeshLib
{
class Mesh;
class NodePartitionedMesh;

/// Function computes all information necessary to transform the MeshLib::Mesh
/// mesh into a NodePartitionedMesh subdomain mesh. Additional to the usual
/// name, nodes and elements, the following information is computed:
/// - a vector of global node ids of the subdomain mesh,
/// - the number of global nodes (the sum of number of nodes of all subdomain
/// meshes)
/// - the number of regular nodes (the sum of number of non-ghost nodes of all
/// subdomain meshes)
/// - a vector containing the number of regular base nodes per rank
/// - a vector containing the number of regular higher order nodes per rank
std::unique_ptr<MeshLib::NodePartitionedMesh>
transformMeshToNodePartitionedMesh(NodePartitionedMesh const* const bulk_mesh,
                                   Mesh const* const subdomain_mesh);
}  // namespace MeshLib
