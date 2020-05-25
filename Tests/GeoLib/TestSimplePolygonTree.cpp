/**
 * @file TestSimplePolygonTree.cpp
 * @author
 * @date Apr 02, 2014
 * @brief Test functionality of class SimplePolygonTree.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    CreatePolygonTreesTest()
    {
        // create points and construct polygon
        pnts_.push_back(new GeoLib::Point(0.0,-1.0,0.0));
        pnts_.push_back(new GeoLib::Point(1.0,0.0,0.0));
        pnts_.push_back(new GeoLib::Point(0.0,1.0,0.0));
        pnts_.push_back(new GeoLib::Point(-1.0,0.0,0.0));

        pnts_.push_back(new GeoLib::Point(-0.9,0.0,0.0));
        pnts_.push_back(new GeoLib::Point(0.75,-0.1,0.0));
        pnts_.push_back(new GeoLib::Point(-0.75,-0.1,0.0));

        // create closed polylines
        GeoLib::Polyline ply0(pnts_);
        ply0.addPoint(3);
        ply0.addPoint(1);
        ply0.addPoint(2);
        ply0.addPoint(3);

        GeoLib::Polyline ply1(pnts_);
        ply1.addPoint(0);
        ply1.addPoint(1);
        ply1.addPoint(3);
        ply1.addPoint(0);

        GeoLib::Polyline ply2(pnts_);
        ply2.addPoint(0);
        ply2.addPoint(1);
        ply2.addPoint(2);
        ply2.addPoint(3);
        ply2.addPoint(0);

        GeoLib::Polyline ply3(pnts_);
        ply3.addPoint(4);
        ply3.addPoint(5);
        ply3.addPoint(6);
        ply3.addPoint(4);

        // create polygons
        p0_ = new GeoLib::Polygon(ply0);
        p1_ = new GeoLib::Polygon(ply1);
        p2_ = new GeoLib::Polygon(ply2);
        p3_ = new GeoLib::Polygon(ply3);
    }

    ~CreatePolygonTreesTest() override
    {
        delete p0_;
        delete p1_;
        delete p2_;
        delete p3_;
        for (auto& pnt_ : pnts_)
        {
            delete pnt_;
        }
    }

protected:
    std::vector<GeoLib::Point*> pnts_;
    GeoLib::Polygon* p0_{nullptr};
    GeoLib::Polygon* p1_{nullptr};
    GeoLib::Polygon* p2_{nullptr};
    GeoLib::Polygon* p3_{nullptr};
};

TEST_F(CreatePolygonTreesTest, P0AndP1)
{
    auto* pt0(new GeoLib::SimplePolygonTree(p0_, nullptr));
    auto* pt1(new GeoLib::SimplePolygonTree(p1_, nullptr));

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
    auto* pt0(new GeoLib::SimplePolygonTree(p0_, nullptr));
    auto* pt1(new GeoLib::SimplePolygonTree(p1_, nullptr));
    auto* pt2(new GeoLib::SimplePolygonTree(p2_, nullptr));

    std::list<GeoLib::SimplePolygonTree*> pt_list;
    pt_list.push_back(pt0);
    pt_list.push_back(pt1);
    pt_list.push_back(pt2);
    ASSERT_EQ(3u, pt_list.size());

    ASSERT_FALSE(p0_->isPolylineInPolygon(*p1_));
    ASSERT_FALSE(p0_->isPolylineInPolygon(*p2_));
    ASSERT_FALSE(p1_->isPolylineInPolygon(*p0_));
    ASSERT_FALSE(p1_->isPolylineInPolygon(*p2_));
    ASSERT_TRUE(p2_->isPolylineInPolygon(*p0_));
    ASSERT_TRUE(p2_->isPolylineInPolygon(*p1_));

    createPolygonTrees(pt_list);
    ASSERT_EQ(1u, pt_list.size());
    ASSERT_EQ(2u, (*(pt_list.begin()))->getNumberOfChildren());
    std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

TEST_F(CreatePolygonTreesTest, P0AndP1AndP2AndP3)
{
    auto* pt0(new GeoLib::SimplePolygonTree(p0_, nullptr));
    auto* pt1(new GeoLib::SimplePolygonTree(p1_, nullptr));
    auto* pt2(new GeoLib::SimplePolygonTree(p2_, nullptr));
    auto* pt3(new GeoLib::SimplePolygonTree(p3_, nullptr));

    std::list<GeoLib::SimplePolygonTree*> pt_list;
    pt_list.push_back(pt0);
    pt_list.push_back(pt1);
    pt_list.push_back(pt2);
    pt_list.push_back(pt3);
    ASSERT_EQ(4u, pt_list.size());

    ASSERT_FALSE(p0_->isPolylineInPolygon(*p1_));
    ASSERT_FALSE(p0_->isPolylineInPolygon(*p2_));
    ASSERT_FALSE(p1_->isPolylineInPolygon(*p0_));
    ASSERT_FALSE(p1_->isPolylineInPolygon(*p2_));
    ASSERT_TRUE(p2_->isPolylineInPolygon(*p0_));
    ASSERT_TRUE(p2_->isPolylineInPolygon(*p1_));
    ASSERT_TRUE(p2_->isPolylineInPolygon(*p3_));
    ASSERT_TRUE(p1_->isPolylineInPolygon(*p3_));

    createPolygonTrees(pt_list);
    ASSERT_EQ(1u, pt_list.size());
    ASSERT_EQ(2u, (*(pt_list.begin()))->getNumberOfChildren());
    std::for_each(pt_list.begin(), pt_list.end(), std::default_delete<GeoLib::SimplePolygonTree>());
}

