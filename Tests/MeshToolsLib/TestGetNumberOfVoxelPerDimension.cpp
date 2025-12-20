// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <array>

#include "MathLib/Point3d.h"
#include "MeshToolsLib/MeshGenerators/VoxelGridFromMesh.h"

TEST(MeshToolsLib, getNumberOfVoxelPerDimension)
{
    std::array<double, 3> const cellsize = {1, 2, 5};
    std::array<double, 3> const ranges = {10, 10, 10};
    std::array<std::size_t, 3> const dims = MeshToolsLib::MeshGenerator::
        VoxelGridFromMesh::getNumberOfVoxelPerDimension(ranges, cellsize);

    ASSERT_EQ(10, dims[0]);
    ASSERT_EQ(5, dims[1]);
    ASSERT_EQ(2, dims[2]);
}

TEST(MeshToolsLib, getNumberOfVoxelPerDimensionMinMax)
{
    std::array<double, 3> const cellsize = {1, 2, 5};
    std::array<double, 3> const ranges = {0, 0, 10};
    std::array<std::size_t, 3> const dims = MeshToolsLib::MeshGenerator::
        VoxelGridFromMesh::getNumberOfVoxelPerDimension(ranges, cellsize);

    ASSERT_EQ(1, dims[0]);
    ASSERT_EQ(1, dims[1]);
    ASSERT_EQ(2, dims[2]);
}

TEST(MeshToolsLib, getNumberOfVoxelPerDimensionCellSize)
{
    std::array<double, 3> const cellsize = {0, 2, 5};
    std::array<double, 3> const ranges = {1, 2, 5};
    ASSERT_ANY_THROW(MeshToolsLib::MeshGenerator::VoxelGridFromMesh::
                         getNumberOfVoxelPerDimension(ranges, cellsize));
}