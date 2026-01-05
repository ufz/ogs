// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
