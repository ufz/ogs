/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "GeoLib/PointVec.h"
#include "GeoLib/Polyline.h"

struct GeoLibCreatePolyline : public testing::Test
{
    GeoLibCreatePolyline()
    {
        points.push_back(new GeoLib::Point{0.0, 0.0, 0.0});
        points.push_back(new GeoLib::Point{1.0, 1.0, 1.0});
        points.push_back(new GeoLib::Point{0.5, 0.5, 0.5});
        std::string point_vec_name("test");
        point_vec = std::make_unique<GeoLib::PointVec>(
            point_vec_name, std::move(points), GeoLib::PointVec::NameIdMap{});
    }

    std::vector<GeoLib::Point*> points;
    std::unique_ptr<GeoLib::PointVec> point_vec;
};

TEST_F(GeoLibCreatePolyline, Empty)
{
    auto polyline = GeoLib::createPolyline(*point_vec, {});
    ASSERT_EQ(0, polyline->getNumberOfPoints());
}

TEST_F(GeoLibCreatePolyline, OnePoint)
{
    auto polyline = GeoLib::createPolyline(*point_vec, {1});
    ASSERT_EQ(1, polyline->getNumberOfPoints());
}

TEST_F(GeoLibCreatePolyline, NonExistingID)
{
    auto polyline = GeoLib::createPolyline(*point_vec, {3});
    ASSERT_EQ(0, polyline->getNumberOfPoints());
}

TEST_F(GeoLibCreatePolyline, ExistingAndNonExistingIDs)
{
    auto polyline = GeoLib::createPolyline(*point_vec, {0, 3, 1, 6, 2});
    ASSERT_EQ(3, polyline->getNumberOfPoints());
    ASSERT_TRUE(polyline->isPointIDInPolyline(0));
    ASSERT_FALSE(polyline->isPointIDInPolyline(3));
    ASSERT_TRUE(polyline->isPointIDInPolyline(1));
    ASSERT_FALSE(polyline->isPointIDInPolyline(6));
    ASSERT_TRUE(polyline->isPointIDInPolyline(2));
}
