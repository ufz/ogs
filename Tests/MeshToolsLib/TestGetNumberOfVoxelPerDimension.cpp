/**
 * \file
 * \date 2023-05-31
 * \brief Tests for getNumberOfVoxelPerDimension()
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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