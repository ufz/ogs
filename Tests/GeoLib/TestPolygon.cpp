/**
 * @file TestPolygon.cpp
 * @author
 * @date Nov 27, 2013
 * @brief Test functionality of class Polygon.
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "Point.h"
#include "Polygon.h"

using namespace GeoLib;

/**
 * Polygon:
 *     3
 *    / \
 *   /   \
 *  2     4
 *  |     |
 *  |     |
 *  1     5
 *   \   /
 *    \ /
 *     0
 */


class IsPntInPolygonTest : public testing::Test
{
public:
	IsPntInPolygonTest() :
		_polygon(nullptr)
	{
		// create points and construct polygon
		_pnts.push_back(new Point(1.0,0.0,0.0));
		_pnts.push_back(new Point(0.0,1.0,0.0));
		_pnts.push_back(new Point(0.0,2.0,0.0));
		_pnts.push_back(new Point(1.0,3.0,0.0));
		_pnts.push_back(new Point(2.0,2.0,0.0));
		_pnts.push_back(new Point(2.0,1.0,0.0));

		// create closed polyline
		Polyline ply(_pnts);
		ply.addPoint(0);
		ply.addPoint(1);
		ply.addPoint(2);
		ply.addPoint(3);
		ply.addPoint(4);
		ply.addPoint(5);
		ply.addPoint(0);

		// create polygon
		_polygon = new Polygon(ply);
	}

	~IsPntInPolygonTest()
	{
		delete _polygon;
		for (std::size_t k(0); k<_pnts.size(); k++)
			delete _pnts[k];
	}

protected:
	std::vector<Point*> _pnts;
	Polygon *_polygon;
};

TEST_F(IsPntInPolygonTest, CheckCorners)
{
	for (std::size_t k(0); k<_pnts.size(); k++)
		EXPECT_TRUE(_polygon->isPntInPolygon(*_pnts[k]));
}

TEST_F(IsPntInPolygonTest, CheckPointsRestOnPolygonEdges)
{
	for (std::size_t k(0); k<_pnts.size()-1; k++) {
		for (double t(0.0); t<1.0; t+=0.001) {
			EXPECT_TRUE(_polygon->isPntInPolygon(
					Point(
							(*_pnts[k])[0] + t*((*_pnts[k+1])[0]-(*_pnts[k])[0]),
							(*_pnts[k])[1] + t*((*_pnts[k+1])[1]-(*_pnts[k])[1]),
							(*_pnts[k])[2] + t*((*_pnts[k+1])[2]-(*_pnts[k])[2])
					))
			);
		}
	}
}

TEST_F(IsPntInPolygonTest, CheckInnerPoints)
{
	ASSERT_TRUE(_polygon->isPntInPolygon(Point(1.0,1.0,0.0)));
	ASSERT_TRUE(_polygon->isPntInPolygon(Point(0.5,1.0,0.0)));
}

TEST_F(IsPntInPolygonTest, CheckOuterPoints)
{
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(1.0-std::numeric_limits<float>::epsilon(),0.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(1.0+std::numeric_limits<float>::epsilon(),0.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(0.0-std::numeric_limits<float>::epsilon(),1.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(0.0-std::numeric_limits<float>::epsilon(),2.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(1.0-std::numeric_limits<float>::epsilon(),3.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(1.0+std::numeric_limits<float>::epsilon(),3.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(2.0+std::numeric_limits<float>::epsilon(),2.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(2.0+std::numeric_limits<float>::epsilon(),1.0,0.0)));
}
