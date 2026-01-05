// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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