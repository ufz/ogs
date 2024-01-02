/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

namespace MeshLib
{
class Mesh;
}

namespace MeshToolsLib::MeshGenerator::AddFaultToVoxelGrid
{
bool addFaultToVoxelGrid(MeshLib::Mesh* mesh,
                         MeshLib::Mesh const* fault,
                         int const fault_id);
// test if input mesh is voxel grid
bool isVoxelGrid(MeshLib::Mesh const& mesh);
}  // namespace MeshToolsLib::MeshGenerator::AddFaultToVoxelGrid