/**
 * \file
 * \author Thomas Fischer
 * \date Jan 7, 2013
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <ctime>
#include <numeric>

#include "GeoLib/Polyline.h"

struct GeoLibPolyline : public testing::Test
{
public:
    GeoLibPolyline() : polyline(points)
    {
        points.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
        points.push_back(new GeoLib::Point(1.0, 0.0, 0.0));
        points.push_back(new GeoLib::Point(0.5, 0.5, 0.0));
        points.push_back(new GeoLib::Point(1.0, 0.5, 0.0));
        points.push_back(new GeoLib::Point(-1.0, 0.0, 0.0));
        points.push_back(new GeoLib::Point(0.0, 0.5, 0.0));
    }

    ~GeoLibPolyline() override
    {
        for (auto& point : points)
        {
            delete point;
        }
    }

protected:
    std::vector<GeoLib::Point*> points;
    GeoLib::Polyline polyline;
};

TEST_F(GeoLibPolyline, Empty)
{
    ASSERT_EQ(std::size_t(0), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_FALSE(polyline.isPointIDInPolyline(0));
}

TEST_F(GeoLibPolyline, TryToAddPointWithToLargeID)
{
    ASSERT_FALSE(polyline.addPoint(points.size()));
    ASSERT_EQ(std::size_t(0), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_FALSE(polyline.isPointIDInPolyline(points.size()));
}

TEST_F(GeoLibPolyline, OnePoint)
{
    ASSERT_TRUE(polyline.addPoint(0));
    // checking properties of polyline with one point
    ASSERT_EQ(std::size_t(1), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(0));
}

TEST_F(GeoLibPolyline, TwoIdenticalPoints)
{
    ASSERT_TRUE(polyline.addPoint(0));
    ASSERT_FALSE(polyline.addPoint(0));
    ASSERT_EQ(std::size_t(1), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(0));
    ASSERT_FALSE(polyline.isPointIDInPolyline(1));
}

TEST_F(GeoLibPolyline, TwoPoints)
{
    ASSERT_TRUE(polyline.addPoint(0));
    ASSERT_TRUE(polyline.addPoint(1));

    // checking properties of polyline with two points
    ASSERT_EQ(std::size_t(2), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(1));
}

TEST_F(GeoLibPolyline, Remove)
{
    ASSERT_TRUE(polyline.addPoint(0));
    ASSERT_TRUE(polyline.addPoint(1));
    ASSERT_TRUE(polyline.addPoint(2));

    // checking properties of polyline with three points
    ASSERT_EQ(std::size_t(3), polyline.getNumberOfPoints());
    ASSERT_TRUE(polyline.isPointIDInPolyline(0));
    ASSERT_TRUE(polyline.isPointIDInPolyline(1));
    ASSERT_TRUE(polyline.isPointIDInPolyline(2));

    // check removing an existing point
    polyline.removePoint(1);
    ASSERT_EQ(std::size_t(2), polyline.getNumberOfPoints());
    ASSERT_TRUE(polyline.isPointIDInPolyline(0));
    ASSERT_FALSE(polyline.isPointIDInPolyline(1));
    ASSERT_TRUE(polyline.isPointIDInPolyline(2));

    // check removing a non-existing point
    polyline.removePoint(polyline.getNumberOfPoints() + 1);
    ASSERT_EQ(std::size_t(2), polyline.getNumberOfPoints());
    ASSERT_TRUE(polyline.isPointIDInPolyline(0));
    ASSERT_FALSE(polyline.isPointIDInPolyline(1));
    ASSERT_TRUE(polyline.isPointIDInPolyline(2));
}

TEST_F(GeoLibPolyline, Insert)
{
    ASSERT_TRUE(polyline.addPoint(0));
    ASSERT_TRUE(polyline.addPoint(2));

    // inserting point in the middle
    ASSERT_TRUE(polyline.insertPoint(1, 1));
    ASSERT_EQ(std::size_t(3), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(2));
    ASSERT_TRUE(polyline.isPointIDInPolyline(0));
    ASSERT_TRUE(polyline.isPointIDInPolyline(1));

    // inserting point at the end
    ASSERT_TRUE(polyline.insertPoint(3, 3));
    ASSERT_EQ(std::size_t(4), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(3));

    // inserting point at the beginning
    ASSERT_TRUE(polyline.insertPoint(0, 4));
    ASSERT_EQ(std::size_t(5), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(4));

    // inserting point in the middle
    ASSERT_TRUE(polyline.insertPoint(2, 5));
    ASSERT_EQ(std::size_t(6), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(5));

    // insert point at non-existing position
    ASSERT_FALSE(polyline.insertPoint(polyline.getNumberOfPoints() + 1, 0));
    ASSERT_EQ(std::size_t(6), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(5));

    // insert non-existing point at existing position
    ASSERT_FALSE(
        polyline.insertPoint(polyline.getNumberOfPoints() - 1, points.size()));
    ASSERT_EQ(std::size_t(6), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(5));

    // insert non-existing point at non-existing position
    ASSERT_FALSE(
        polyline.insertPoint(polyline.getNumberOfPoints() + 1, points.size()));
    ASSERT_EQ(std::size_t(6), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    ASSERT_TRUE(polyline.isPointIDInPolyline(5));
}

TEST_F(GeoLibPolyline, Close)
{
    std::vector<std::size_t> point_ids(points.size());
    std::iota(point_ids.begin(), point_ids.end(), 0);

    for (auto point_id : point_ids)
    {
        ASSERT_TRUE(polyline.addPoint(point_id));
    }

    // close polyline; this adds the first polyline point to the end
    polyline.closePolyline();
    ASSERT_EQ(points.size() + 1, polyline.getNumberOfPoints());
    ASSERT_TRUE(polyline.isClosed());

    // remove last point -> polyline is not closed!
    polyline.removePoint(polyline.getNumberOfPoints() - 1);
    ASSERT_EQ(point_ids.size(), polyline.getNumberOfPoints());
    ASSERT_FALSE(polyline.isClosed());
    polyline.closePolyline();
    ASSERT_TRUE(polyline.isClosed());
}

TEST_F(GeoLibPolyline, CountSegments)
{
    polyline.addPoint(0);
    polyline.addPoint(1);
    polyline.addPoint(2);

    std::size_t segment_cnt(0);
    for (auto seg_it(polyline.begin()); seg_it != polyline.end(); ++seg_it)
    {
        ++segment_cnt;
    }
    ASSERT_EQ(polyline.getNumberOfSegments(), segment_cnt);
}
