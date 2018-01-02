/**
 * @file TestPolyline.cpp
 * @author Thomas Fischer
 * @date Jan 7, 2013
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>

#include "gtest/gtest.h"

#include "GeoLib/Polyline.h"

TEST(GeoLib, PolylineTest)
{
    std::vector<GeoLib::Point*> ply_pnts;
    GeoLib::Polyline ply(ply_pnts);

    // checking properties of empty polyline
    ASSERT_EQ(std::size_t(0), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_FALSE(ply.isPointIDInPolyline(0));
    ASSERT_EQ(0.0, ply.getLength(0));

    ply_pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    ply.addPoint(0);
    // checking properties of polyline with one point
    ASSERT_EQ(std::size_t(1), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_TRUE(ply.isPointIDInPolyline(0));
    ASSERT_EQ(0.0, ply.getLength(0));

    ply_pnts.push_back(new GeoLib::Point(1.0, 0.0, 0.0));
    ply.addPoint(1);
    // checking properties of polyline with two points
    ASSERT_EQ(std::size_t(2), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_TRUE(ply.isPointIDInPolyline(1));
    ASSERT_EQ(1.0, ply.getLength(1));

    // checking properties of polyline with two points
    ply_pnts.push_back(new GeoLib::Point(0.5, 0.5, 0.0));
    ply.addPoint(2);
    ASSERT_EQ(std::size_t(3), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_TRUE(ply.isPointIDInPolyline(2));
    ASSERT_TRUE(fabs(ply.getLength(2) - (1.0 + sqrt(0.5))) < std::numeric_limits<double>::epsilon());

    // checking remove
    ply.removePoint(1);
    ASSERT_EQ(std::size_t(2), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_FALSE(ply.isPointIDInPolyline(1));
    ASSERT_TRUE(fabs(ply.getLength(1) - sqrt(0.5)) < std::numeric_limits<double>::epsilon());

    // inserting point in the middle
    ply.insertPoint(1,1);
    ASSERT_EQ(std::size_t(3), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_TRUE(ply.isPointIDInPolyline(2));
    ASSERT_TRUE(fabs(ply.getLength(2) - (1.0 + sqrt(0.5))) < std::numeric_limits<double>::epsilon());

    // inserting point at the end
    ply_pnts.push_back(new GeoLib::Point(1.0, 0.5, 0.0));
    ply.insertPoint(3,3);
    ASSERT_EQ(std::size_t(4), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_TRUE(ply.isPointIDInPolyline(3));
    ASSERT_TRUE(fabs(ply.getLength(3) - (1.0 + sqrt(0.5) + 0.5)) < std::numeric_limits<double>::epsilon());

    // inserting point at the beginning
    ply_pnts.push_back(new GeoLib::Point(-1.0, 0.0, 0.0));
    ply.insertPoint(0,4);
    ASSERT_EQ(std::size_t(5), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_TRUE(ply.isPointIDInPolyline(4));
    ASSERT_TRUE(fabs(ply.getLength(4) - (1.0 + 1.0 + sqrt(0.5) + 0.5)) <
                std::numeric_limits<double>::epsilon());

    // inserting point in the middle
    ply_pnts.push_back(new GeoLib::Point(0.0, 0.5, 0.0));
    ply.insertPoint(2,5);
    ASSERT_EQ(std::size_t(6), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ASSERT_TRUE(ply.isPointIDInPolyline(5));
    ASSERT_TRUE(fabs(ply.getLength(5) - (1.0 + 0.5 + sqrt(1.25) + sqrt(0.5) + 0.5)) < std::numeric_limits<double>::epsilon());

    // close polyline
    ply.closePolyline();
    ASSERT_EQ(std::size_t(7), ply.getNumberOfPoints());
    ASSERT_TRUE(ply.isClosed());

    // remove last point -> polyline is not closed!
    ply.removePoint(6);
    ASSERT_EQ(std::size_t(6), ply.getNumberOfPoints());
    ASSERT_FALSE(ply.isClosed());
    ply.closePolyline();
    ASSERT_TRUE(ply.isClosed());

    std::size_t segment_cnt(0);
    for (auto seg_it(ply.begin()); seg_it != ply.end(); ++seg_it)
    {
        ++segment_cnt;
    }
    ASSERT_EQ(ply.getNumberOfSegments(), segment_cnt);

    for (auto & ply_pnt : ply_pnts)
        delete ply_pnt;
}
