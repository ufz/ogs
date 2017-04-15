/**
 * @file TestSimplePolygonTree.cpp
 * @author
 * @date Apr 02, 2014
 * @brief Test functionality of class SimplePolygonTree.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>

#include "gtest/gtest.h"

#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/SimplePolygonTree.h"

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
        _pnts.push_back(new GeoLib::Point(0.0,-1.0,0.0));
        _pnts.push_back(new GeoLib::Point(1.0,0.0,0.0));
        _pnts.push_back(new GeoLib::Point(0.0,1.0,0.0));
        _pnts.push_back(new GeoLib::Point(-1.0,0.0,0.0));

        _pnts.push_back(new GeoLib::Point(-0.9,0.0,0.0));
        _pnts.push_back(new GeoLib::Point(0.75,-0.1,0.0));
        _pnts.push_back(new GeoLib::Point(-0.75,-0.1,0.0));

        // create closed polylines
        GeoLib::Polyline ply0(_pnts);
        ply0.addPoint(3);
        ply0.addPoint(1);
        ply0.addPoint(2);
        ply0.addPoint(3);

        GeoLib::Polyline ply1(_pnts);
        ply1.addPoint(0);
        ply1.addPoint(1);
        ply1.addPoint(3);
        ply1.addPoint(0);

        GeoLib::Polyline ply2(_pnts);
        ply2.addPoint(0);
        ply2.addPoint(1);
        ply2.addPoint(2);
        ply2.addPoint(3);
        ply2.addPoint(0);

        GeoLib::Polyline ply3(_pnts);
        ply3.addPoint(4);
        ply3.addPoint(5);
        ply3.addPoint(6);
        ply3.addPoint(4);

        // create polygons
        _p0 = new GeoLib::Polygon(ply0);
        _p1 = new GeoLib::Polygon(ply1);
        _p2 = new GeoLib::Polygon(ply2);
        _p3 = new GeoLib::Polygon(ply3);
    }

    ~CreatePolygonTreesTest()
    {
        delete _p0;
        delete _p1;
        delete _p2;
        delete _p3;
        for (auto & _pnt : _pnts)
            delete _pnt;
    }

protected:
    std::vector<GeoLib::Point*> _pnts;
    GeoLib::Polygon *_p0;
    GeoLib::Polygon *_p1;
    GeoLib::Polygon *_p2;
    GeoLib::Polygon *_p3;
};

TEST_F(CreatePolygonTreesTest, P0AndP1)
{
    auto* pt0(new GeoLib::SimplePolygonTree(_p0, nullptr));
    auto* pt1(new GeoLib::SimplePolygonTree(_p1, nullptr));

    std::list<GeoLib::SimplePolygonTree*> pt_list;
    pt_list.push_back(pt0);
    pt_list.push_back(pt1);
    ASSERT_EQ(2u, pt_list.size());

    createPolygonTrees(pt_list);
    ASSERT_EQ(2u, pt_list.size());
    std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

TEST_F(CreatePolygonTreesTest, P0AndP1AndP2)
{
    auto* pt0(new GeoLib::SimplePolygonTree(_p0, nullptr));
    auto* pt1(new GeoLib::SimplePolygonTree(_p1, nullptr));
    auto* pt2(new GeoLib::SimplePolygonTree(_p2, nullptr));

    std::list<GeoLib::SimplePolygonTree*> pt_list;
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
    ASSERT_EQ(2u, (*(pt_list.begin()))->getNumberOfChildren());
    std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

TEST_F(CreatePolygonTreesTest, P0AndP1AndP2AndP3)
{
    auto* pt0(new GeoLib::SimplePolygonTree(_p0, nullptr));
    auto* pt1(new GeoLib::SimplePolygonTree(_p1, nullptr));
    auto* pt2(new GeoLib::SimplePolygonTree(_p2, nullptr));
    auto* pt3(new GeoLib::SimplePolygonTree(_p3, nullptr));

    std::list<GeoLib::SimplePolygonTree*> pt_list;
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
    ASSERT_EQ(2u, (*(pt_list.begin()))->getNumberOfChildren());
    std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

