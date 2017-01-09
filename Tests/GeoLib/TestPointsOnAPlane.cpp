/**
 * \file   TestPointsOnAPlane.cpp
 * \author Karsten Rink
 * \date   2014-02-26
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <ctime>
#include <random>

#include <gtest/gtest.h>

#include "MathLib/GeometricBasics.h"
#include "GeoLib/Point.h"

void testAllPossibilities(
    GeoLib::Point const& a,
    GeoLib::Point const& b,
    GeoLib::Point const& c,
    GeoLib::Point const& d,
    bool expected)
{
    ASSERT_EQ(expected, MathLib::isCoplanar(a, b, c, d));
    ASSERT_EQ(expected, MathLib::isCoplanar(a, b, d, c));
    ASSERT_EQ(expected, MathLib::isCoplanar(a, c, b, d));
    ASSERT_EQ(expected, MathLib::isCoplanar(a, c, d, b));
    ASSERT_EQ(expected, MathLib::isCoplanar(a, d, b, c));
    ASSERT_EQ(expected, MathLib::isCoplanar(a, d, c, b));

    ASSERT_EQ(expected, MathLib::isCoplanar(b, a, c, d));
    ASSERT_EQ(expected, MathLib::isCoplanar(b, a, d, c));
    ASSERT_EQ(expected, MathLib::isCoplanar(b, c, a, d));
    ASSERT_EQ(expected, MathLib::isCoplanar(b, c, d, a));
    ASSERT_EQ(expected, MathLib::isCoplanar(b, d, a, c));
    ASSERT_EQ(expected, MathLib::isCoplanar(b, d, c, a));

    ASSERT_EQ(expected, MathLib::isCoplanar(c, b, a, d));
    ASSERT_EQ(expected, MathLib::isCoplanar(c, b, d, a));
    ASSERT_EQ(expected, MathLib::isCoplanar(c, a, b, d));
    ASSERT_EQ(expected, MathLib::isCoplanar(c, a, d, b));
    ASSERT_EQ(expected, MathLib::isCoplanar(c, d, b, a));
    ASSERT_EQ(expected, MathLib::isCoplanar(c, d, a, b));

    ASSERT_EQ(expected, MathLib::isCoplanar(d, b, c, a));
    ASSERT_EQ(expected, MathLib::isCoplanar(d, b, a, c));
    ASSERT_EQ(expected, MathLib::isCoplanar(d, c, b, a));
    ASSERT_EQ(expected, MathLib::isCoplanar(d, c, a, b));
    ASSERT_EQ(expected, MathLib::isCoplanar(d, a, b, c));
    ASSERT_EQ(expected, MathLib::isCoplanar(d, a, c, b));
}

TEST(GeoLib, TestPointsOnAPlane)
{
    // 2d case
    GeoLib::Point a(0,0,0);
    GeoLib::Point b(1,0,0);
    GeoLib::Point c(0,1,0);
    GeoLib::Point d(1,1,0);
    testAllPossibilities(a, b, c, d, true);

    // disturbe z coordinate of point d
    d[2] = 1e6*std::numeric_limits<double>::epsilon();
    testAllPossibilities(a, b, c, d, true);

    // disturbe z coordinate of point d
    d[2] = 1e-9;
    testAllPossibilities(a, b, c, d, true);
    d[2] = 1e-6;
    testAllPossibilities(a, b, c, d, true);
    d[2] = 1e-5;
    testAllPossibilities(a, b, c, d, false);

    // irregular case: a = d = ORIGIN
    d = GeoLib::Point(0.0, 0.0, 0.0);
    testAllPossibilities(a, b, c, d, true);

    // irregular case: a = b = d = ORIGIN
    b = GeoLib::Point(0.0, 0.0, 0.0);
    testAllPossibilities(a, b, c, d, true);

    // irregular case: a = b = c = d
    c = GeoLib::Point(0.0, 0.0, 0.0);
    testAllPossibilities(a, b, c, d, true);

    // 2d case with fixed z
    a = GeoLib::Point(0.0, 0.0, 0.3);
    b = GeoLib::Point(1e-5, 0.0, 0.3);
    c = GeoLib::Point(0.0, 1e-5, 0.3);
    d = GeoLib::Point(1e-5, 1e-5, 0.3);
    testAllPossibilities(a, b, c, d, true);

    // a,b,c with random coordinates,
    std::uniform_real_distribution<double> distri(-1e10, 1e10);
    std::default_random_engine re;
    for (std::size_t k(0); k<1000; k++) {
        a = GeoLib::Point(distri(re), distri(re), distri(re));
        b = GeoLib::Point(distri(re), distri(re), distri(re));
        c = GeoLib::Point(distri(re), distri(re), distri(re));
        // d such that it is in the plane
        d = GeoLib::Point(b[0]+c[0]-a[0], b[1]+c[1]-a[1], b[2]+c[2]-a[2]);
        testAllPossibilities(a, b, c, d, true);
        // d such that it is not in the plane
        d[2] = std::numeric_limits<double>::epsilon() +
            (1 + std::numeric_limits<double>::epsilon())*(b[2]+c[2]);
        testAllPossibilities(a, b, c, d, false);
    }
}

