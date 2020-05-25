/**
 * \file
 * \date Nov 27, 2013
 * \brief Test functionality of class Polygon.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include "GeoLib/Point.h"
#include "GeoLib/LineSegment.h"
#include "GeoLib/Polygon.h"

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
    PolygonTest()
    {
        // create points and construct polygon
        pnts_.push_back(new GeoLib::Point( 0.0, 0.0,0.0)); // 0
        pnts_.push_back(new GeoLib::Point(-2.0, 2.0,0.0)); // 1
        pnts_.push_back(new GeoLib::Point(-2.0, 4.0,0.0)); // 2
        pnts_.push_back(new GeoLib::Point(-1.0, 2.0,0.0)); // 3
        pnts_.push_back(new GeoLib::Point( 0.0, 4.0,0.0)); // 4
        pnts_.push_back(new GeoLib::Point( 1.0, 2.0,0.0)); // 5
        pnts_.push_back(new GeoLib::Point( 2.0, 4.0,0.0)); // 6
        pnts_.push_back(new GeoLib::Point( 2.0, 2.0,0.0)); // 7

        // create closed polyline
        GeoLib::Polyline ply(pnts_);
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
        polygon_ = new GeoLib::Polygon(ply);
    }

    ~PolygonTest() override
    {
        delete polygon_;
        for (auto& pnt_ : pnts_)
        {
            delete pnt_;
        }
    }

protected:
    std::vector<GeoLib::Point*> pnts_;
    GeoLib::Polygon* polygon_{nullptr};
};

TEST_F(PolygonTest, isPntInPolygonCheckCorners)
{
    for (auto& pnt_ : pnts_)
    {
        EXPECT_TRUE(polygon_->isPntInPolygon(*pnt_));
    }
}

TEST_F(PolygonTest, isPntInPolygonCheckPointsRestOnPolygonEdges)
{
    for (std::size_t k(0); k<pnts_.size()-1; k++) {
        double t = 0;
        while (t < 1.0)
        {
            EXPECT_TRUE(polygon_->isPntInPolygon(GeoLib::Point(
                (*pnts_[k])[0] + t * ((*pnts_[k + 1])[0] - (*pnts_[k])[0]),
                (*pnts_[k])[1] + t * ((*pnts_[k + 1])[1] - (*pnts_[k])[1]),
                (*pnts_[k])[2] + t * ((*pnts_[k + 1])[2] - (*pnts_[k])[2]))));
            t += 0.001;
        }
    }
}

TEST_F(PolygonTest, isPntInPolygonCheckInnerPoints)
{
    ASSERT_TRUE(polygon_->isPntInPolygon(GeoLib::Point(1.0,1.0,0.0)));
    ASSERT_TRUE(polygon_->isPntInPolygon(GeoLib::Point(0.5,1.0,0.0)));
}

TEST_F(PolygonTest, isPntInPolygonCheckOuterPoints)
{
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        0.0-std::numeric_limits<float>::epsilon(),0.0,0.0)));
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        -2.0-std::numeric_limits<float>::epsilon(),2.0,0.0)));
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        -2.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        -1.0, 2.0+std::numeric_limits<float>::epsilon(),0.0)));
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        0.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        1.0,2.0+std::numeric_limits<float>::epsilon(),0.0)));
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        2.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
    ASSERT_FALSE(polygon_->isPntInPolygon(GeoLib::Point(
        2.0+std::numeric_limits<float>::epsilon(),2.0,0.0)));
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
TEST_F(PolygonTest, containsSegment)
{
    { // test segment (2,6)
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(polygon_->getPoint(2)),
            const_cast<GeoLib::Point*>(polygon_->getPoint(6))};
        ASSERT_FALSE(polygon_->containsSegment(segment));
    }

    // test all segments of polygon
    for (auto && segment_it : *polygon_)
    {
        EXPECT_TRUE(polygon_->containsSegment(segment_it));
    }

    { // 70
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(polygon_->getPoint(7)),
            const_cast<GeoLib::Point*>(polygon_->getPoint(0))};
        ASSERT_TRUE(polygon_->containsSegment(segment));
    }

    { // test segment (3,5)
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(polygon_->getPoint(3)),
            const_cast<GeoLib::Point*>(polygon_->getPoint(5))};
        ASSERT_TRUE(polygon_->containsSegment(segment));
    }

    { // test segment (1,7)
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(polygon_->getPoint(1)),
            const_cast<GeoLib::Point*>(polygon_->getPoint(7))};
        ASSERT_TRUE(polygon_->containsSegment(segment));
    }

    { // test segment (1,4)
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(polygon_->getPoint(1)),
            const_cast<GeoLib::Point*>(polygon_->getPoint(4))};
        ASSERT_FALSE(polygon_->containsSegment(segment));
    }
}

TEST_F(PolygonTest, isPolylineInPolygon)
{
    // create a test polyline
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(-2.0,4.0,0.0)); // 2
    pnts.push_back(new GeoLib::Point( 2.0,4.0,0.0)); // 6
    GeoLib::Polyline outer_ply(pnts);
    outer_ply.addPoint(0);
    outer_ply.addPoint(1);
    ASSERT_FALSE(polygon_->isPolylineInPolygon(outer_ply));
    for (auto& pnt : pnts)
    {
        delete pnt;
    }
    pnts.clear();

    pnts.push_back(new GeoLib::Point(-1.0,2.0,0.0)); // 3
    pnts.push_back(new GeoLib::Point( 1.0,2.0,0.0)); // 5
    GeoLib::Polyline inner_ply(pnts);
    inner_ply.addPoint(0);
    inner_ply.addPoint(1);
    ASSERT_TRUE(polygon_->isPolylineInPolygon(inner_ply));
    for (auto& pnt : pnts)
    {
        delete pnt;
    }
}

TEST_F(PolygonTest, CopyConstructor)
{
    GeoLib::Polygon polygon_copy(*polygon_);
    ASSERT_EQ(polygon_->getNumberOfSegments(), polygon_copy.getNumberOfSegments());

    // Check if all line segments of the original polygon are contained in the
    // copy
    for (auto const& segment : *polygon_)
    {
        ASSERT_TRUE(polygon_copy.containsSegment(segment));
    }
}
