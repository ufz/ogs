/**
 * @file TestIsPointInTriangle.cpp
 * @date 2014-05-27
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "MathLib/Point3d.h"
#include "MathLib/GeometricBasics.h"

TEST(MathLib, IsPointInTriangle)
{
    MathLib::Point3d const a(std::array<double, 3>{{0.0, 0.0, 0.0}});
    MathLib::Point3d const b(std::array<double, 3>{{100.0, 0.0, 0.0}});
    MathLib::Point3d const c(std::array<double, 3>{{0.0, 100.0, 0.0}});

    // check point on corner points of triangle
    MathLib::Point3d q(a);
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = b;
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = c;
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c));

    // check points on edges of triangle
    q = MathLib::Point3d(std::array<double, 3>{
        {0.5 * (b[0] + a[0]), 0.5 * (b[1] + a[1]), 0.5 * (b[2] + a[2])}});
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = MathLib::Point3d(std::array<double, 3>{
        {0.5 * (c[0] + a[0]), 0.5 * (c[1] + a[1]), 0.5 * (c[2] + a[2])}});
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = MathLib::Point3d(std::array<double, 3>{
        {0.5 * (c[0] + b[0]), 0.5 * (c[1] + b[1]), 0.5 * (c[2] + b[2])}});
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c));

    // check points inside
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 0.0}});
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 1e-10}});
    EXPECT_TRUE(MathLib::gaussPointInTriangle(q, a, b, c));
    // here is a higher eps value needed for the second algorithm
    EXPECT_TRUE(MathLib::barycentricPointInTriangle(q, a, b, c, 1e-5));

    // check points outside
    q = MathLib::Point3d(std::array<double, 3>{{-0.1, 0.1, 0.0}});
    EXPECT_FALSE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_FALSE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 0.0005}});
    EXPECT_FALSE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_FALSE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 0.0001}});
    EXPECT_FALSE(MathLib::gaussPointInTriangle(
        q, a, b, c, std::numeric_limits<double>::epsilon()));
    EXPECT_FALSE(MathLib::barycentricPointInTriangle(
        q, a, b, c, std::numeric_limits<double>::epsilon()));
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 0.000001}});
    EXPECT_FALSE(MathLib::gaussPointInTriangle(
        q, a, b, c, std::numeric_limits<double>::epsilon()));
    EXPECT_FALSE(MathLib::barycentricPointInTriangle(
        q, a, b, c, std::numeric_limits<double>::epsilon()));
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 1e-7}});
    EXPECT_FALSE(MathLib::gaussPointInTriangle(
        q, a, b, c, std::numeric_limits<double>::epsilon()));
    EXPECT_FALSE(MathLib::barycentricPointInTriangle(
        q, a, b, c, std::numeric_limits<double>::epsilon()));
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 0.001}});
    EXPECT_FALSE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_FALSE(MathLib::barycentricPointInTriangle(q, a, b, c));
    q = MathLib::Point3d(std::array<double, 3>{{0.1, 0.1, 0.1}});
    EXPECT_FALSE(MathLib::gaussPointInTriangle(q, a, b, c));
    EXPECT_FALSE(MathLib::barycentricPointInTriangle(q, a, b, c));
}

