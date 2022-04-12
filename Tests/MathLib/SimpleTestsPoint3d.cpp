/**
 * @file TestPoint3d.cpp
 * @date Mar 13, 2014
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <array>

#include "MathLib/Point3d.h"

using namespace MathLib;

TEST(MathLib, Point3dConstructors)
{
    // *** test default constructors for different DIM values
    Point3d p0;
    // test coordinates of default constructed point
    ASSERT_EQ(0.0, p0[0]);
    ASSERT_EQ(0.0, p0[1]);
    ASSERT_EQ(0.0, p0[2]);

    // *** test copy constructor
    Point3d p0_copy(p0);
    // test equality of coordinates
    ASSERT_EQ(p0[0], p0_copy[0]);
    ASSERT_EQ(p0[1], p0_copy[1]);
    ASSERT_EQ(p0[2], p0_copy[2]);

    // *** test constructor taking std::array
    std::array array = {0.0, 1.0, 2.0};
    Point3d p3(array);
    ASSERT_EQ(0.0, p3[0]);
    ASSERT_EQ(1.0, p3[1]);
    ASSERT_EQ(2.0, p3[2]);
}

TEST(MathLib, Point3dOperators)
{
    Point3d p;
    // access operator
    p[0] = 1.0;
    p[1] = 3.0;
    p[2] = 5.0;
    ASSERT_EQ(1.0, p[0]);
    ASSERT_EQ(3.0, p[1]);
    ASSERT_EQ(5.0, p[2]);
}
