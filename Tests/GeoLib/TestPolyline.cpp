/**
 * @file TestPolyline.cpp
 * @author Thomas Fischer
 * @date Jan 7, 2013
 *
 * @copyright
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "Polyline.h"

using namespace GeoLib;

TEST(GeoLib, PolylineConstructorTest)
{
	std::vector<Point*> ply_pnts;
	Polyline ply(ply_pnts);

	// checking properties of empty polyline
	ASSERT_EQ(ply.getNumberOfPoints(), 0);
	ASSERT_FALSE(ply.isClosed());
	ASSERT_FALSE(ply.isPointIDInPolyline(0));
	ASSERT_EQ(ply.getLength(0), 0.0);

	ply_pnts.push_back(new Point(0.0, 0.0, 0.0));
	ply.addPoint(0);
	// checking properties of polyline with one point
	ASSERT_EQ(ply.getNumberOfPoints(), 1);
	ASSERT_FALSE(ply.isClosed());
	ASSERT_TRUE(ply.isPointIDInPolyline(0));
	ASSERT_EQ(ply.getLength(0), 0.0);

	ply_pnts.push_back(new Point(1.0, 0.0, 0.0));
	ply.addPoint(1);
	// checking properties of polyline with two points
	ASSERT_EQ(ply.getNumberOfPoints(), 2);
	ASSERT_FALSE(ply.isClosed());
	ASSERT_TRUE(ply.isPointIDInPolyline(1));
	ASSERT_EQ(ply.getLength(1), 1);

	ply_pnts.push_back(new Point(0.5, 0.5, 0.0));
	ply.addPoint(2);
	ASSERT_EQ(ply.getNumberOfPoints(), 3);
	ASSERT_FALSE(ply.isClosed());
	ASSERT_TRUE(ply.isPointIDInPolyline(2));
	ASSERT_TRUE(fabs(ply.getLength(2) - (1.0 + sqrt(0.5))) < std::numeric_limits<double>::min());
}
