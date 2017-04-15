/**
 * @file TestPolygon.cpp
 * @author
 * @date Nov 27, 2013
 * @brief Test functionality of class Polygon.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
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
    PolygonTest() :
        _polygon(nullptr)
    {
        // create points and construct polygon
        _pnts.push_back(new GeoLib::Point( 0.0, 0.0,0.0)); // 0
        _pnts.push_back(new GeoLib::Point(-2.0, 2.0,0.0)); // 1
        _pnts.push_back(new GeoLib::Point(-2.0, 4.0,0.0)); // 2
        _pnts.push_back(new GeoLib::Point(-1.0, 2.0,0.0)); // 3
        _pnts.push_back(new GeoLib::Point( 0.0, 4.0,0.0)); // 4
        _pnts.push_back(new GeoLib::Point( 1.0, 2.0,0.0)); // 5
        _pnts.push_back(new GeoLib::Point( 2.0, 4.0,0.0)); // 6
        _pnts.push_back(new GeoLib::Point( 2.0, 2.0,0.0)); // 7

        // create closed polyline
        GeoLib::Polyline ply(_pnts);
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
        _polygon = new GeoLib::Polygon(ply);
    }

    ~PolygonTest() override
    {
        delete _polygon;
        for (auto & _pnt : _pnts)
            delete _pnt;
    }

protected:
    std::vector<GeoLib::Point*> _pnts;
    GeoLib::Polygon *_polygon;
};

TEST_F(PolygonTest, isPntInPolygonCheckCorners)
{
    for (auto & _pnt : _pnts)
        EXPECT_TRUE(_polygon->isPntInPolygon(*_pnt));
}

TEST_F(PolygonTest, isPntInPolygonCheckPointsRestOnPolygonEdges)
{
    for (std::size_t k(0); k<_pnts.size()-1; k++) {
        for (double t(0.0); t<1.0; t+=0.001) {
            EXPECT_TRUE(_polygon->isPntInPolygon(
                GeoLib::Point(
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
    ASSERT_TRUE(_polygon->isPntInPolygon(GeoLib::Point(1.0,1.0,0.0)));
    ASSERT_TRUE(_polygon->isPntInPolygon(GeoLib::Point(0.5,1.0,0.0)));
}

TEST_F(PolygonTest, isPntInPolygonCheckOuterPoints)
{
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
        0.0-std::numeric_limits<float>::epsilon(),0.0,0.0)));
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
        -2.0-std::numeric_limits<float>::epsilon(),2.0,0.0)));
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
        -2.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
        -1.0, 2.0+std::numeric_limits<float>::epsilon(),0.0)));
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
        0.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
        1.0,2.0+std::numeric_limits<float>::epsilon(),0.0)));
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
        2.0-std::numeric_limits<float>::epsilon(),4.0,0.0)));
    ASSERT_FALSE(_polygon->isPntInPolygon(GeoLib::Point(
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
            const_cast<GeoLib::Point*>(_polygon->getPoint(2)),
            const_cast<GeoLib::Point*>(_polygon->getPoint(6))};
        ASSERT_FALSE(_polygon->containsSegment(segment));
    }

    // test all segments of polygon
    for (auto && segment_it : *_polygon)
    {
        EXPECT_TRUE(_polygon->containsSegment(segment_it));
    }

    { // 70
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(_polygon->getPoint(7)),
            const_cast<GeoLib::Point*>(_polygon->getPoint(0))};
        ASSERT_TRUE(_polygon->containsSegment(segment));
    }

    { // test segment (3,5)
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(_polygon->getPoint(3)),
            const_cast<GeoLib::Point*>(_polygon->getPoint(5))};
        ASSERT_TRUE(_polygon->containsSegment(segment));
    }

    { // test segment (1,7)
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(_polygon->getPoint(1)),
            const_cast<GeoLib::Point*>(_polygon->getPoint(7))};
        ASSERT_TRUE(_polygon->containsSegment(segment));
    }

    { // test segment (1,4)
        GeoLib::LineSegment const segment{
            const_cast<GeoLib::Point*>(_polygon->getPoint(1)),
            const_cast<GeoLib::Point*>(_polygon->getPoint(4))};
        ASSERT_FALSE(_polygon->containsSegment(segment));
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
    ASSERT_FALSE(_polygon->isPolylineInPolygon(outer_ply));
    for (auto & pnt : pnts)
        delete pnt;
    pnts.clear();

    pnts.push_back(new GeoLib::Point(-1.0,2.0,0.0)); // 3
    pnts.push_back(new GeoLib::Point( 1.0,2.0,0.0)); // 5
    GeoLib::Polyline inner_ply(pnts);
    inner_ply.addPoint(0);
    inner_ply.addPoint(1);
    ASSERT_TRUE(_polygon->isPolylineInPolygon(inner_ply));
    for (auto & pnt : pnts)
        delete pnt;
}

TEST_F(PolygonTest, CopyConstructor)
{
    GeoLib::Polygon polygon_copy(*_polygon);
    ASSERT_EQ(_polygon->getNumberOfSegments(), polygon_copy.getNumberOfSegments());

    // Check if all line segments of the original polygon are contained in the
    // copy
    for (auto const& segment : *_polygon)
        ASSERT_TRUE(polygon_copy.containsSegment(segment));
}
