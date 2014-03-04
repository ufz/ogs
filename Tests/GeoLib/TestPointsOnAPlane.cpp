/**
 * \file   TestPointsOnAPlane.cpp
 * \author Karsten Rink
 * \date   2014-02-26
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "gtest/gtest.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/Point.h"

TEST(GeoLib, TestPointsOnAPlane)
{
	GeoLib::Point a(0,0,0);
	GeoLib::Point b(1,0,0);
	GeoLib::Point c(0,1,0);
	GeoLib::Point d(1,1,0);

	bool result = GeoLib::pointsOnAPlane(a, b, c, d);
	ASSERT_EQ(result, true);

	d = GeoLib::Point(1,1,0.1);
	result = GeoLib::dividedByPlane(a, b, c, d);
	ASSERT_EQ(result, false);
}
