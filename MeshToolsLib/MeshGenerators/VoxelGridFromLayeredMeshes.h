// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace MeshLib
{
class Mesh;
}

namespace MathLib
{
class Point3d;
}

namespace MeshToolsLib::MeshGenerators::VoxelFromLayeredMeshes
{

/**
 * Constructs a VoxelGrid with a specified resolution of a list of
 * layered meshes.
 * \param extent  The axis-aligned boundary-box (AABB) which is initially
 * defined.
 * \param layers                 Containing all the meshes that will be used
 * to create the VoxelGrid.
 * \param cellsize       Contains the resolution of the Voxel (i.e.
 * length in x-,y-,z-directions)
 * \param dilate  A flag to set dilate. If dilate is True all Voxels which
 * are not fully covered by the meshes are included in the resulting
 * VoxelGrid.
 */
std::unique_ptr<MeshLib::Mesh> createVoxelFromLayeredMesh(
    std::pair<MathLib::Point3d, MathLib::Point3d>& extent,
    std::vector<MeshLib::Mesh const*> const& layers,
    std::array<double, 3> const cellsize,
    bool const dilate);
}  // namespace MeshToolsLib::MeshGenerators::VoxelFromLayeredMeshes
