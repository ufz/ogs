/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <numeric>

#include "CreateTestPoints.h"
#include "GeoLib/PointVec.h"
#include "GeoLib/Polyline.h"

struct GeoLibResetPolylinePointIDs : public testing::Test
{
    GeoLibResetPolylinePointIDs()
    {
        auto const limits =
            std::array<double, 6>{-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
        auto random_points = createRandomPoints(51, limits);

        point_vec =
            std::make_unique<GeoLib::PointVec>("test_points",
                                               std::move(random_points),
                                               GeoLib::PointVec::NameIdMap{});
    }

    std::unique_ptr<GeoLib::PointVec> point_vec;
};

TEST_F(GeoLibResetPolylinePointIDs, Empty)
{
    auto polyline = GeoLib::createPolyline(*point_vec, {});
    std::vector<std::size_t> empty_mapping{};
    ASSERT_THROW(GeoLib::resetPointIDs(*(polyline.get()), empty_mapping),
                 std::runtime_error);
}

TEST_F(GeoLibResetPolylinePointIDs, AllPoints)
{
    // create polyline using all points
    std::vector<std::size_t> ids(point_vec->size());
    std::iota(ids.begin(), ids.end(), 0ul);
    auto polyline = GeoLib::createPolyline(*point_vec, std::move(ids));

    // create mapping that doesn't change anything
    std::vector<std::size_t> identity_mapping(point_vec->size());
    std::iota(identity_mapping.begin(), identity_mapping.end(), 0ul);

    GeoLib::resetPointIDs(*(polyline.get()), identity_mapping);

    for (std::size_t k = 0; k < polyline->getNumberOfPoints(); ++k)
    {
        ASSERT_EQ(k, polyline->getPointID(k));
    }
}

TEST_F(GeoLibResetPolylinePointIDs, ReverseAllPoints)
{
    // create polyline using all points
    std::vector<std::size_t> ids(point_vec->size());
    std::iota(ids.begin(), ids.end(), 0ul);
    auto polyline = GeoLib::createPolyline(*point_vec, std::move(ids));

    // create mapping that reverses the sequence of points
    std::vector<std::size_t> reverse_mapping(point_vec->size());
    std::generate(reverse_mapping.begin(), reverse_mapping.end(),
                  [n = point_vec->size()]() mutable { return --n; });

    GeoLib::resetPointIDs(*(polyline.get()), reverse_mapping);

    auto const n = polyline->getNumberOfPoints();
    for (std::size_t k = 0; k < n; ++k)
    {
        ASSERT_EQ(n - 1 - k, polyline->getPointID(k));
    }
}

TEST_F(GeoLibResetPolylinePointIDs, HalfNumberOfPoints)
{
    // create polyline using half number of points (from the end of PointsVec)
    std::vector<std::size_t> ids(point_vec->size() / 2);
    std::iota(ids.begin(), ids.end(), point_vec->size() / 2);
    auto polyline = GeoLib::createPolyline(*point_vec, std::move(ids));

    // create mapping
    std::vector<std::size_t> mapping(point_vec->size());
    std::iota(mapping.begin(),
              mapping.begin() + mapping.size() / 2,
              mapping.size() / 2);
    std::iota(mapping.begin() + mapping.size() / 2, mapping.end(), 0ul);

    GeoLib::resetPointIDs(*(polyline.get()), mapping);

    auto const n = polyline->getNumberOfPoints();
    for (std::size_t k = 0; k < n; ++k)
    {
        ASSERT_EQ(k, polyline->getPointID(k));
    }
}
