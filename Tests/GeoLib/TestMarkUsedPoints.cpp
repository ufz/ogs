/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <numeric>

#include "BaseLib/Algorithm.h"
#include "CreateTestPoints.h"
#include "GeoLib/PointVec.h"
#include "GeoLib/Polyline.h"

struct GeoLibMarkUsedPoints : public testing::Test
{
    GeoLibMarkUsedPoints()
    {
        auto const limits =
            std::array<double, 6>{-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
        auto random_points = createRandomPoints(51, limits);

        point_vec = std::make_unique<GeoLib::PointVec>(
            "test_points", std::move(random_points),
            GeoLib::PointVec::NameIdMap{});
    }

    std::unique_ptr<GeoLib::PointVec> point_vec;
};

TEST_F(GeoLibMarkUsedPoints, InvalidInputVector)
{
    auto polyline = GeoLib::createPolyline(*point_vec, {});
    std::vector<bool> empty_used_points{};
    ASSERT_THROW(GeoLib::markUsedPoints(*(polyline.get()), empty_used_points),
                 std::runtime_error);
}

TEST_F(GeoLibMarkUsedPoints, AllPoints)
{
    // create polyline using all points
    std::vector<std::size_t> ids(point_vec->size());
    std::iota(ids.begin(), ids.end(), 0ul);
    auto polyline = GeoLib::createPolyline(*point_vec, std::move(ids));

    std::vector<bool> used_points(point_vec->size(), false);

    GeoLib::markUsedPoints(*(polyline.get()), used_points);

    ASSERT_TRUE(std::all_of(used_points.begin(),
                            used_points.end(),
                            [](auto const bool_value) { return bool_value; }));
}

TEST_F(GeoLibMarkUsedPoints, HalfNumberOfPoints)
{
    // create polyline using half number of points (from the end of PointsVec)
    std::vector<std::size_t> ids(point_vec->size() / 2);
    std::iota(ids.begin(), ids.end(), point_vec->size() / 2);
    auto polyline = GeoLib::createPolyline(*point_vec, std::move(ids));

    std::vector<bool> used_points(point_vec->size(), false);

    GeoLib::markUsedPoints(*(polyline.get()), used_points);

    ASSERT_EQ(polyline->getNumberOfPoints(),
              std::count(used_points.begin(), used_points.end(), true));
}
