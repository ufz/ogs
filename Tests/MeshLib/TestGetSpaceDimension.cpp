/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/Algorithm.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/GetSpaceDimension.h"

// nodes on x axis
TEST(MeshLib, GetSpaceDimensionOnXAxis)
{
    std::vector<MeshLib::Node*> nodes(2);
    nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
    nodes[1] = new MeshLib::Node(0.2, 0.0, 0.0, 1);

    auto const expected_space_dimension = 1u;
    ASSERT_EQ(expected_space_dimension, MeshLib::getSpaceDimension(nodes));
    BaseLib::cleanupVectorElements(nodes);
}

// nodes on y axis
TEST(MeshLib, GetSpaceDimensionOnYAxis)
{
    std::vector<MeshLib::Node*> nodes(2);
    nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
    nodes[1] = new MeshLib::Node(0.0, 2.0, 0.0, 1);

    auto const expected_space_dimension = 2u;
    ASSERT_EQ(expected_space_dimension, MeshLib::getSpaceDimension(nodes));
    BaseLib::cleanupVectorElements(nodes);
}
// nodes on z axis
TEST(MeshLib, GetSpaceDimensionOnZAxis)
{
    std::vector<MeshLib::Node*> nodes(2);
    nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0, 0);
    nodes[1] = new MeshLib::Node(0.0, 0.0, 2.0, 1);

    auto const expected_space_dimension = 3u;
    ASSERT_EQ(expected_space_dimension, MeshLib::getSpaceDimension(nodes));
    BaseLib::cleanupVectorElements(nodes);
}

// nodes on x-y plane
TEST(MeshLib, GetSpaceDimensionOnXYPlane)
{
    std::vector<MeshLib::Node*> nodes(2);
    nodes[0] = new MeshLib::Node(-1.0, 0.0, 0.0, 0);
    nodes[1] = new MeshLib::Node(0.0, 2.0, 0.0, 1);

    auto const expected_space_dimension = 2u;
    ASSERT_EQ(expected_space_dimension, MeshLib::getSpaceDimension(nodes));
    BaseLib::cleanupVectorElements(nodes);
}

// nodes on y-z plane
TEST(MeshLib, GetSpaceDimensionOnYZPlane)
{
    std::vector<MeshLib::Node*> nodes(2);
    nodes[0] = new MeshLib::Node(0.0, -1.0, .0, 0);
    nodes[1] = new MeshLib::Node(0.0, 0.0, 2.0, 1);

    auto const expected_space_dimension = 3u;
    ASSERT_EQ(expected_space_dimension, MeshLib::getSpaceDimension(nodes));
    BaseLib::cleanupVectorElements(nodes);
}

// nodes on x-z plane
TEST(MeshLib, GetSpaceDimensionOnXZPlane)
{
    std::vector<MeshLib::Node*> nodes(2);
    nodes[0] = new MeshLib::Node(-1.0, 0.0, .0, 0);
    nodes[1] = new MeshLib::Node(0.0, 0.0, 2.0, 1);

    auto const expected_space_dimension = 3u;
    ASSERT_EQ(expected_space_dimension, MeshLib::getSpaceDimension(nodes));
    BaseLib::cleanupVectorElements(nodes);
}

// nodes in 3D space
TEST(MeshLib, GetSpaceDimensionIn3DSpace)
{
    std::vector<MeshLib::Node*> nodes(3);
    nodes[0] = new MeshLib::Node(-1.0, 1.0, .0, 0);
    nodes[1] = new MeshLib::Node(0.0, 3.0, 2.0, 1);
    nodes[2] = new MeshLib::Node(1.0, 0.0, 2.0, 1);

    auto const expected_space_dimension = 3u;
    ASSERT_EQ(expected_space_dimension, MeshLib::getSpaceDimension(nodes));
    BaseLib::cleanupVectorElements(nodes);
}
