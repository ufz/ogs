/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"

using namespace MeshLib;

TEST(MeshLib, UniqueMeshId)
{
    // Create first mesh and get the current counter value.
    Mesh m0("first", std::vector<Node*>(), std::vector<Element*>());
    std::size_t const counter_value = m0.getID();

    EXPECT_GE(counter_value, 0u);

    //
    // Test mesh counter increments.
    //
    Mesh* m1 = new Mesh("second", std::vector<Node*>(), std::vector<Element*>());
    ASSERT_EQ(counter_value + std::size_t(1), m1->getID());

    Mesh m2("third", std::vector<Node*>(), std::vector<Element*>());
    ASSERT_EQ(counter_value + std::size_t(2), m2.getID());

    delete m1;
    ASSERT_EQ(counter_value + std::size_t(2), m2.getID());

    Mesh m3("fourth", std::vector<Node*>(), std::vector<Element*>());
    ASSERT_EQ(counter_value + std::size_t(3), m3.getID());

    // Copy mesh keeps also increments the counter.
    Mesh m4 = m0;
    ASSERT_EQ(counter_value + std::size_t(4), m4.getID());

}

