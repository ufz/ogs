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
 *  2     4     6
 *  |\   / \   /|
 *  | \ /   \ / |
 *  1  3     5  7
 *   \         /
 *    \       /
 *     \     /
 *      \   /
 *       \ /
 *        0
 */


class PolygonTest : public testing::Test
{
public:
	PolygonTest() :
		_polygon(nullptr)
	{
		// create points and construct polygon
		_pnts.push_back(new Point( 0.0, 0.0,0.0)); // 0
		_pnts.push_back(new Point(-2.0, 2.0,0.0)); // 1
		_pnts.push_back(new Point(-2.0, 4.0,0.0)); // 2
		_pnts.push_back(new Point(-1.0, 2.0,0.0)); // 3
		_pnts.push_back(new Point( 0.0, 4.0,0.0)); // 4
		_pnts.push_back(new Point( 1.0, 2.0,0.0)); // 5
		_pnts.push_back(new Point( 2.0, 4.0,0.0)); // 6
		_pnts.push_back(new Point( 2.0, 2.0,0.0)); // 7

		// create closed polyline
		Polyline ply(_pnts);
		ply.addPoint(0);
		ply.addPoint(1);
		ply.addPoint(2);
		ply.addPoint(3);
		ply.addPoint(4);
		ply.addPoint(5);
		ply.addPoint(6);
		ply.addPoint(7);
		ply.addPoint(0);

		// create polygon
		_polygon = new Polygon(ply);
	}

	~PolygonTest()
	{
		delete _polygon;
		for (std::size_t k(0); k<_pnts.size(); k++)
			delete _pnts[k];
	}

protected:
	std::vector<Point*> _pnts;
	Polygon *_polygon;
};

TEST_F(PolygonTest, isPntInPolygonCheckCorners)
{
	for (std::size_t k(0); k<_pnts.size(); k++)
		EXPECT_TRUE(_polygon->isPntInPolygon(*_pnts[k]));
}

TEST_F(PolygonTest, isPntInPolygonCheckPointsRestOnPolygonEdges)
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

TEST_F(PolygonTest, isPntInPolygonCheckInnerPoints)
{
	ASSERT_TRUE(_polygon->isPntInPolygon(Point(1.0,1.0,0.0)));
	ASSERT_TRUE(_polygon->isPntInPolygon(Point(0.5,1.0,0.0)));
}

TEST_F(PolygonTest, isPntInPolygonCheckOuterPoints)
{
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(0.0-std::numeric_limits<float>::epsilon(),0.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(-2.0-std::numeric_limits<float>::epsilon(),2.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(-2.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(-1.0, 2.0+std::numeric_limits<float>::epsilon(),0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(0.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(1.0,2.0+std::numeric_limits<float>::epsilon(),0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(2.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
	ASSERT_FALSE(_polygon->isPntInPolygon(Point(2.0+std::numeric_limits<float>::epsilon(),2.0,0.0)));
}

/**
 *  2     4     6
 *  |\   / \   /|
 *  | \ /   \ / |
 *  1  3     5  7
 *   \         /
 *    \       /
 *     \     /
 *      \   /
 *       \ /
 *        0
 * 0 = (0,0), 1=(-2,2), 2=(-2,4), 3=(-1,2), 4=(0,4), 5=(1,2), 6=(2,4), 7=(2,2)
 */

TEST_F(PolygonTest, isPolylineInPolygon)
{
	// create a test polyline
	std::vector<GeoLib::Point*> pnts;
	pnts.push_back(new GeoLib::Point(-2.0,4.0,0.0)); // 2
	pnts.push_back(new GeoLib::Point( 2.0,4.0,0.0)); // 6
	GeoLib::Polyline outer_ply(pnts);
	outer_ply.addPoint(0);
	outer_ply.addPoint(1);
	ASSERT_FALSE(_polygon->isPolylineInPolygon(outer_ply));
	for (std::size_t k(0); k<pnts.size(); k++)
		delete pnts[k];
	pnts.clear();

	pnts.push_back(new GeoLib::Point(-1.0,2.0,0.0)); // 3
	pnts.push_back(new GeoLib::Point( 1.0,2.0,0.0)); // 5
	GeoLib::Polyline inner_ply(pnts);
	inner_ply.addPoint(0);
	inner_ply.addPoint(1);
	ASSERT_TRUE(_polygon->isPolylineInPolygon(inner_ply));
	for (std::size_t k(0); k<pnts.size(); k++)
		delete pnts[k];
}
