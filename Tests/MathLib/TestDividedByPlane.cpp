/**
 * \file   TestDividedByPlane.cpp
 * \author Karsten Rink
 * \date   2014-02-26
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>

#include "MathLib/GeometricBasics.h"
#include "MathLib/Point3d.h"

TEST(MathLib, TestDividedByPlane)
{
    // xy plane
    MathLib::Point3d a(std::array<double,3>{{0, 0, 0}});
    MathLib::Point3d b(std::array<double,3>{{1, 0, 0}});
    MathLib::Point3d c(std::array<double,3>{{0, 1, 0}});
    MathLib::Point3d d(std::array<double,3>{{1, 1, 0}});

    bool result = MathLib::dividedByPlane(a, d, b, c);
    ASSERT_TRUE(result);

    result = MathLib::dividedByPlane(b, c, a, d);
    ASSERT_TRUE(result);

    d = MathLib::Point3d(std::array<double,3>{{0.1, 0.1, 0}});
    result = MathLib::dividedByPlane(b, c, a, d);
    ASSERT_FALSE(result);

    // xz plane
    c = MathLib::Point3d(std::array<double,3>{{0, 0, 1}});
    d = MathLib::Point3d(std::array<double,3>{{1, 0, 1}});
    result = MathLib::dividedByPlane(a, d, b, c);
    ASSERT_TRUE(result);

    // yz plane
    b = MathLib::Point3d(std::array<double,3>{{0, 1, 0}});
    d = MathLib::Point3d(std::array<double,3>{{0, 1, 1}});
    result = MathLib::dividedByPlane(a, d, b, c);
    ASSERT_TRUE(result);
}
