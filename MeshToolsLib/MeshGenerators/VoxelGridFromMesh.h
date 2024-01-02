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

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <array>
#include <memory>

namespace MeshLib
{
class Mesh;
}

namespace MathLib
{
class Point3d;
}
namespace MeshToolsLib::MeshGenerator::VoxelGridFromMesh
{
static std::string const cell_id_name = "CellIds";

/// getNumberOfVoxelPerDimension is used to calculate how many voxel fit into
/// a bounding box. For this calculation the difference of min and max point of
/// the bounding box is divided by the cell size, for every dimension. The
/// calculation is restricted to work only with positive values for the cell
/// size. If the difference between min and max is zero, we assign one voxel for
/// the respective dimension.
std::array<std::size_t, 3> getNumberOfVoxelPerDimension(
    std::array<double, 3> const& ranges, std::array<double, 3> const& cellsize);

std::vector<int> assignCellIds(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                               MathLib::Point3d const& min,
                               std::array<std::size_t, 3> const& dims,
                               std::array<double, 3> const& cellsize);

// grid has to contain a PropertyVector with the name 'CellIds'
bool removeUnusedGridCells(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                           std::unique_ptr<MeshLib::Mesh>& grid);
// map the cell data of mesh to voxelgrid
void mapMeshArraysOntoGrid(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                           std::unique_ptr<MeshLib::Mesh> const& grid);
}  // namespace MeshToolsLib::MeshGenerator::VoxelGridFromMesh
