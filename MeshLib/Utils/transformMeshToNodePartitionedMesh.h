/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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

// compute all information necessary to transform the MeL::Mesh mesh into a
// NodePartitionedMesh subdomain mesh
std::unique_ptr<MeshLib::NodePartitionedMesh>
transformMeshToNodePartitionedMesh(NodePartitionedMesh const* const bulk_mesh,
                                   Mesh const* const subdomain_mesh);
}  // namespace MeshLib
