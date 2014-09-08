/**
 * @file TestSimplePolygonTree.cpp
 * @author
 * @date Apr 02, 2014
 * @brief Test functionality of class SimplePolygonTree.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>

#include "gtest/gtest.h"

#include "Point.h"
#include "Polygon.h"
#include "SimplePolygonTree.h"

using namespace GeoLib;

/**
 *       2
 *      / \
 *     /   \
 *    /     \
 *   /       \
 *  3---------1
 *   \  6--5 /
 *    \ \ / /
 *     \ 4 /
 *      \ /
 *       0
 *
 * Polygons: P0: 3,1,2,3 / P1: 0,1,3,0 / P2: 0,1,2,3,0
 *           P3: 4,5,6,4
 */


class CreatePolygonTreesTest : public testing::Test
{
public:
	CreatePolygonTreesTest() :
		_p0(nullptr), _p1(nullptr), _p2(nullptr), _p3(nullptr)
	{
		// create points and construct polygon
		_pnts.push_back(new Point(0.0,-1.0,0.0));
		_pnts.push_back(new Point(1.0,0.0,0.0));
		_pnts.push_back(new Point(0.0,1.0,0.0));
		_pnts.push_back(new Point(-1.0,0.0,0.0));

		_pnts.push_back(new Point(-0.9,0.0,0.0));
		_pnts.push_back(new Point(0.75,-0.1,0.0));
		_pnts.push_back(new Point(-0.75,-0.1,0.0));

		// create closed polylines
		Polyline ply0(_pnts);
		ply0.addPoint(3);
		ply0.addPoint(1);
		ply0.addPoint(2);
		ply0.addPoint(3);

		Polyline ply1(_pnts);
		ply1.addPoint(0);
		ply1.addPoint(1);
		ply1.addPoint(3);
		ply1.addPoint(0);

		Polyline ply2(_pnts);
		ply2.addPoint(0);
		ply2.addPoint(1);
		ply2.addPoint(2);
		ply2.addPoint(3);
		ply2.addPoint(0);

		Polyline ply3(_pnts);
		ply3.addPoint(4);
		ply3.addPoint(5);
		ply3.addPoint(6);
		ply3.addPoint(4);

		// create polygons
		_p0 = new Polygon(ply0);
		_p1 = new Polygon(ply1);
		_p2 = new Polygon(ply2);
		_p3 = new Polygon(ply3);
	}

	~CreatePolygonTreesTest()
	{
		delete _p0;
		delete _p1;
		delete _p2;
		delete _p3;
		for (std::size_t k(0); k<_pnts.size(); k++)
			delete _pnts[k];
	}

protected:
	std::vector<Point*> _pnts;
	Polygon *_p0;
	Polygon *_p1;
	Polygon *_p2;
	Polygon *_p3;
};

TEST_F(CreatePolygonTreesTest, P0AndP1)
{
	SimplePolygonTree *pt0(new SimplePolygonTree(_p0, nullptr));
	SimplePolygonTree *pt1(new SimplePolygonTree(_p1, nullptr));

	std::list<SimplePolygonTree*> pt_list;
	pt_list.push_back(pt0);
	pt_list.push_back(pt1);
	ASSERT_EQ(2u, pt_list.size());

	createPolygonTrees(pt_list);
	ASSERT_EQ(2u, pt_list.size());
	std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

TEST_F(CreatePolygonTreesTest, P0AndP1AndP2)
{
	SimplePolygonTree *pt0(new SimplePolygonTree(_p0, nullptr));
	SimplePolygonTree *pt1(new SimplePolygonTree(_p1, nullptr));
	SimplePolygonTree *pt2(new SimplePolygonTree(_p2, nullptr));

	std::list<SimplePolygonTree*> pt_list;
	pt_list.push_back(pt0);
	pt_list.push_back(pt1);
	pt_list.push_back(pt2);
	ASSERT_EQ(3u, pt_list.size());

	ASSERT_FALSE(_p0->isPolylineInPolygon(*_p1));
	ASSERT_FALSE(_p0->isPolylineInPolygon(*_p2));
	ASSERT_FALSE(_p1->isPolylineInPolygon(*_p0));
	ASSERT_FALSE(_p1->isPolylineInPolygon(*_p2));
	ASSERT_TRUE(_p2->isPolylineInPolygon(*_p0));
	ASSERT_TRUE(_p2->isPolylineInPolygon(*_p1));

	createPolygonTrees(pt_list);
	ASSERT_EQ(1u, pt_list.size());
	ASSERT_EQ(2u, (*(pt_list.begin()))->getNChilds());
	std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

TEST_F(CreatePolygonTreesTest, P0AndP1AndP2AndP3)
{
	SimplePolygonTree *pt0(new SimplePolygonTree(_p0, nullptr));
	SimplePolygonTree *pt1(new SimplePolygonTree(_p1, nullptr));
	SimplePolygonTree *pt2(new SimplePolygonTree(_p2, nullptr));
	SimplePolygonTree *pt3(new SimplePolygonTree(_p3, nullptr));

	std::list<SimplePolygonTree*> pt_list;
	pt_list.push_back(pt0);
	pt_list.push_back(pt1);
	pt_list.push_back(pt2);
	pt_list.push_back(pt3);
	ASSERT_EQ(4u, pt_list.size());

	ASSERT_FALSE(_p0->isPolylineInPolygon(*_p1));
	ASSERT_FALSE(_p0->isPolylineInPolygon(*_p2));
	ASSERT_FALSE(_p1->isPolylineInPolygon(*_p0));
	ASSERT_FALSE(_p1->isPolylineInPolygon(*_p2));
	ASSERT_TRUE(_p2->isPolylineInPolygon(*_p0));
	ASSERT_TRUE(_p2->isPolylineInPolygon(*_p1));
	ASSERT_TRUE(_p2->isPolylineInPolygon(*_p3));
	ASSERT_TRUE(_p1->isPolylineInPolygon(*_p3));

	createPolygonTrees(pt_list);
	ASSERT_EQ(1u, pt_list.size());
	ASSERT_EQ(2u, (*(pt_list.begin()))->getNChilds());
	std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

