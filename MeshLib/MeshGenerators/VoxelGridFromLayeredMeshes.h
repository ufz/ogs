/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/ProjectPointOnMesh.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"

namespace MeshLib::MeshGenerators::VoxelFromLayeredMeshes
{

/**
 * Constructs a VoxelGrid with a specified resolution of a list of
 * layered meshes.
 * @param extent  The axis-aligned boundary-box (AABB) which is initially
 * defined.
 * @param layers                 Containing all the meshes that will be used
 * to create the VoxelGrid.
 * @param cellsize       Contains the resolution of the Voxel (i.e.
 * length in x-,y-,z-directions)
 * @param dilate  A flag to set dilate. If dilate is True all Voxels which
 * are not fully covered by the meshes are included in the resulting
 * VoxelGrid.
 */
std::unique_ptr<MeshLib::Mesh> createVoxelFromLayeredMesh(
    std::pair<MathLib::Point3d, MathLib::Point3d>& extent,
    std::vector<MeshLib::Mesh const*> const& layers,
    std::array<double, 3> const cellsize,
    bool const dilate);
}  // namespace MeshLib::MeshGenerators::VoxelFromLayeredMeshes