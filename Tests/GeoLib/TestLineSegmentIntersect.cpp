/**
 * @file TestLineSegmentIntersect.cpp
 * @date Dec 2, 2013
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include <random>

#include "gtest/gtest.h"

#include "GeoLib/Point.h"
#include "GeoLib/AnalyticalGeometry.h"

TEST(GeoLib, TestXAxisParallelLineSegmentIntersection2d)
{
    // two intersecting parallel line segments
    GeoLib::Point a(0.0, 0.0, 0.0);
    GeoLib::Point b(1.0, 0.0, 0.0);
    GeoLib::Point c(0.5, 0.0, 0.0);
    GeoLib::Point d(1.5, 0.0, 0.0);
    GeoLib::Point s;

    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestParallelLineSegmentIntersection2d)
{
    // two intersecting parallel line segments
    GeoLib::Point a(0.0, 0.0, 0.0);
    GeoLib::Point b(1.0, 1.0, 0.0);
    GeoLib::Point c(0.5, 0.5, 0.0);
    GeoLib::Point d(1.5, 1.5, 0.0);
    GeoLib::Point s;

    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestNonIntersectingParallelLineSegments2dCaseI)
{
    // two non intersecting parallel line segments
    GeoLib::Point a(0.0, 0.0, 0.0);
    GeoLib::Point b(1.0, 1.0, 0.0);
    GeoLib::Point c(1.5, 1.5, 0.0);
    GeoLib::Point d(2.5, 2.5, 0.0);
    GeoLib::Point s;

    EXPECT_FALSE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestNonIntersectingParallelLineSegments2dCaseII)
{
    // two non intersecting parallel line segments
    GeoLib::Point a(0.0, 1.0, 0.0);
    GeoLib::Point b(1.0, 0.0, 0.0);
    GeoLib::Point c(0.0, -1.0, 0.0);
    GeoLib::Point d(-1.0, 0.0, 0.0);
    GeoLib::Point s;

    EXPECT_FALSE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestIntersectingLineSegments2d)
{
    // two intersecting line segments
    GeoLib::Point a(0.0, 0.0, 0.0);
    GeoLib::Point b(1.0, 1.0, 0.0);
    GeoLib::Point c(0.0, 1.0, 0.0);
    GeoLib::Point d(1.0, 0.0, 0.0);
    GeoLib::Point s;

    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestIntersectingLineSegments3d)
{
    // two intersecting line segments
    GeoLib::Point a(0.0, 0.0, 0.0);
    GeoLib::Point b(1.0, 1.0, 1.0);
    GeoLib::Point c(0.0, 1.0, 0.0);
    GeoLib::Point d(1.0, 0.0, 1.0);
    GeoLib::Point s;

    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));

    // disturb one point a little bit
    double const eps(std::numeric_limits<float>::epsilon());
    d[2] += eps;
    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
    d[2] = 1.0+1e-9;
    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
    d[2] = 1.0+1e-8;
    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
    d[2] = 1.0 + 5e-6;
    EXPECT_FALSE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestParallelLineSegmentIntersection3d)
{
    // two parallel line segments
    GeoLib::Point a(0.0, 0.0, 0.0);
    GeoLib::Point b(1.0, 1.0, 1.0);
    GeoLib::Point c(0.5, 0.5, 0.5);
    GeoLib::Point d(1.5, 1.5, 1.5);
    GeoLib::Point s;

    EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestNonIntersectingParallelLineSegments3d)
{
    // two non intersecting parallel line segments
    GeoLib::Point a(0.0, 0.0, 0.0);
    GeoLib::Point b(1.0, 1.0, 1.0);
    GeoLib::Point c(1.5, 1.5, 1.5);
    GeoLib::Point d(2.5, 2.5, 2.5);
    GeoLib::Point s;
    EXPECT_FALSE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));

    // two axis parallel line segments that do not intersect
    b[0] = 3.0;
    b[1] = 0.0;
    b[2] = 0.0;
    c[0] = 4.0;
    c[1] = 0.0;
    c[2] = 0.0;
    d[0] = 10.0;
    d[1] = 0.0;
    d[2] = 0.0;
    EXPECT_FALSE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));

    // two axis parallel line segments that do not intersect
    a[0] = 1000.0;
    a[1] = 0.0;
    a[2] = 100.0;
    b[0] = 400.0;
    b[1] = 0.0;
    b[2] = 100.0;
    c[0] = 300.0;
    c[1] = 0.0;
    c[2] = 100.0;
    d[0] = 0.0;
    d[1] = 0.0;
    d[2] = 100.0;
    EXPECT_FALSE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
}

TEST(GeoLib, TestRandomLineSegments3d)
{
    std::uniform_real_distribution<double> distribution(-1e10, 1e10);
    std::default_random_engine re;

    // line segments intersect each other in the middle
    for (std::size_t k(0); k<1000; k++) {
        // two intersecting line segments (a,b) and (c,d)
        GeoLib::Point a(distribution(re), distribution(re), distribution(re));
        GeoLib::Point b(distribution(re), distribution(re), distribution(re));
        // direction of (c,d)
        GeoLib::Point w(distribution(re), distribution(re), distribution(re));
        // construct c and d such that (a,b) intersects with (c,d)
        GeoLib::Point c(0.5*(b[0]+a[0]) + w[0],
                0.5*(b[1]+a[1]) + w[1],
                0.5*(b[2]+a[2]) + w[2]);
        GeoLib::Point d(0.5*(b[0]+a[0]) - w[0],
                0.5*(b[1]+a[1]) - w[1],
                0.5*(b[2]+a[2]) - w[2]);
        GeoLib::Point s;

        EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
    }

    // line segment (c,d) intersects (a,b) at end point a
    for (std::size_t k(0); k<1000; k++) {
        // two intersecting line segments (a,b) and (c,d)
        GeoLib::Point a(distribution(re), distribution(re), distribution(re));
        GeoLib::Point b(distribution(re), distribution(re), distribution(re));
        // direction of (c,d)
        GeoLib::Point w(distribution(re), distribution(re), distribution(re));
        // construct c and d such that (a,b) intersects with (c,d)
        GeoLib::Point c(a[0] + w[0], a[1] + w[1], a[2] + w[2]);
        GeoLib::Point d(a[0] - w[0], a[1] - w[1], a[2] - w[2]);
        GeoLib::Point s;

        EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
        EXPECT_TRUE(GeoLib::lineSegmentIntersect({&c,&d}, {&a,&b}, s));
    }

    // line segment (c,d) intersects (a,b) at end point b
    for (std::size_t k(0); k<1000; k++) {
        // two intersecting line segments (a,b) and (c,d)
        GeoLib::Point a(distribution(re), distribution(re), distribution(re));
        GeoLib::Point b(distribution(re), distribution(re), distribution(re));
        // direction of (c,d)
        GeoLib::Point w(distribution(re), distribution(re), distribution(re));
        // construct c and d such that (a,b) intersects with (c,d)
        GeoLib::Point c(b[0] + w[0], b[1] + w[1], b[2] + w[2]);
        GeoLib::Point d(b[0] - w[0], b[1] - w[1], b[2] - w[2]);
        GeoLib::Point s;

        EXPECT_TRUE(GeoLib::lineSegmentIntersect({&a,&b}, {&c,&d}, s));
        EXPECT_TRUE(GeoLib::lineSegmentIntersect({&c,&d}, {&a,&b}, s));
    }
}

