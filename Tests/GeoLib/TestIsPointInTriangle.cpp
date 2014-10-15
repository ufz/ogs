/**
 * @file TestIsPointInTriangle.cpp
 * @date 2014-05-27
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include "AnalyticalGeometry.h"

TEST(GeoLib, IsPointInTriangle)
{
	GeoLib::Point const a(0.0, 0.0, 0.0);
	GeoLib::Point const b(100.0, 0.0, 0.0);
	GeoLib::Point const c(0.0, 100.0, 0.0);

	double const default_eps (std::numeric_limits<double>::epsilon());

	// check point on corner points of triangle
	GeoLib::Point q(a);
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = b;
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = c;
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c));

	// check points on edges of triangle
	q = GeoLib::Point(0.5*(b[0]+a[0]), 0.5*(b[1]+a[1]), 0.5*(b[2]+a[2]));
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = GeoLib::Point(0.5*(c[0]+a[0]), 0.5*(c[1]+a[1]), 0.5*(c[2]+a[2]));
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = GeoLib::Point(0.5*(c[0]+b[0]), 0.5*(c[1]+b[1]), 0.5*(c[2]+b[2]));
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c));

	// check points inside
	q = GeoLib::Point (0.1, 0.1, 0.0);
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = GeoLib::Point (0.1, 0.1, 1e-10);
	EXPECT_TRUE(GeoLib::gaussPointInTriangle(q, a, b, c));
	// here is a higher eps value needed for the second algorithm
	EXPECT_TRUE(GeoLib::barycentricPointInTriangle(q, a, b, c, 1e-5));

	// check points outside
	q = GeoLib::Point (-0.1, 0.1, 0.0);
	EXPECT_FALSE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_FALSE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = GeoLib::Point (0.1, 0.1, 0.0005);
	EXPECT_FALSE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_FALSE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = GeoLib::Point (0.1, 0.1, 0.0001);
	EXPECT_FALSE(GeoLib::gaussPointInTriangle(q, a, b, c, std::numeric_limits<double>::epsilon()));
	EXPECT_FALSE(GeoLib::barycentricPointInTriangle(q, a, b, c, std::numeric_limits<double>::epsilon()));
	q = GeoLib::Point (0.1, 0.1, 0.000001);
	EXPECT_FALSE(GeoLib::gaussPointInTriangle(q, a, b, c, std::numeric_limits<double>::epsilon()));
	EXPECT_FALSE(GeoLib::barycentricPointInTriangle(q, a, b, c, std::numeric_limits<double>::epsilon()));
	q = GeoLib::Point (0.1, 0.1, 1e-7);
	EXPECT_FALSE(GeoLib::gaussPointInTriangle(q, a, b, c, std::numeric_limits<double>::epsilon()));
	EXPECT_FALSE(GeoLib::barycentricPointInTriangle(q, a, b, c, std::numeric_limits<double>::epsilon()));
	q = GeoLib::Point (0.1, 0.1, 0.001);
	EXPECT_FALSE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_FALSE(GeoLib::barycentricPointInTriangle(q, a, b, c));
	q = GeoLib::Point (0.1, 0.1, 0.1);
	EXPECT_FALSE(GeoLib::gaussPointInTriangle(q, a, b, c));
	EXPECT_FALSE(GeoLib::barycentricPointInTriangle(q, a, b, c));
}

