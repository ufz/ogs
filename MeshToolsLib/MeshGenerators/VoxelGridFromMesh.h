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

std::array<std::size_t, 3> getNumberOfVoxelPerDimension(
    std::array<double, 3> const& ranges, std::array<double, 3> const& cellsize);

std::vector<int> assignCellIds(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                               MathLib::Point3d const& min,
                               std::array<std::size_t, 3> const& dims,
                               std::array<double, 3> const& cellsize);

bool removeUnusedGridCells(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                           std::unique_ptr<MeshLib::Mesh>& grid);

void mapMeshArraysOntoGrid(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                           std::unique_ptr<MeshLib::Mesh>& grid);
};  // namespace MeshToolsLib::MeshGenerator::VoxelGridFromMesh