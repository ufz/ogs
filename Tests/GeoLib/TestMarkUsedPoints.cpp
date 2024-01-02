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

#include "BaseLib/Algorithm.h"
#include "CreateTestPoints.h"
#include "GeoLib/PointVec.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

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

struct GeoLibMarkUsedPointsInSurface : public testing::Test
{
    GeoLibMarkUsedPointsInSurface()
    {
        auto const distances = std::array<double, 3>{1.0, 1.0, 1.0};
        number_of_subdivisions = std::array<unsigned, 3>{9, 10, 0};
        point_vec = createPointsInGridArrangement(number_of_subdivisions,
                                                  distances, MathLib::ORIGIN);
        used_points = std::vector<bool>(point_vec.size(), false);
    }

    ~GeoLibMarkUsedPointsInSurface()
    {
        BaseLib::cleanupVectorElements(point_vec);
    }

    std::vector<bool> used_points;
    std::array<unsigned, 3> number_of_subdivisions;
    std::vector<GeoLib::Point*> point_vec;
};

TEST_F(GeoLibMarkUsedPointsInSurface, EmptySurface)
{
    // create empty surface
    auto surface = GeoLib::Surface(point_vec);

    GeoLib::markUsedPoints(surface, used_points);

    ASSERT_EQ(surface.getNumberOfTriangles(),
              std::count(used_points.begin(), used_points.end(), true));
}

TEST_F(GeoLibMarkUsedPointsInSurface, OneTriangle)
{
    // create surface base on one triangle
    auto surface = GeoLib::Surface(point_vec);
    surface.addTriangle(0, 1, number_of_subdivisions[0] + 1);

    GeoLib::markUsedPoints(surface, used_points);

    ASSERT_EQ(3, std::count(used_points.begin(), used_points.end(), true));
}

TEST_F(GeoLibMarkUsedPointsInSurface, TwoTriangles)
{
    // create surface base on two triangles that share two points
    auto surface = GeoLib::Surface(point_vec);
    surface.addTriangle(0, 1, number_of_subdivisions[0] + 1);
    surface.addTriangle(1, number_of_subdivisions[0] + 2,
                        number_of_subdivisions[0] + 1);

    GeoLib::markUsedPoints(surface, used_points);

    ASSERT_EQ(4, std::count(used_points.begin(), used_points.end(), true));
}

TEST_F(GeoLibMarkUsedPointsInSurface, AllPointsUsedInTriangles)
{
    // create surface covering all points
    auto surface = GeoLib::Surface(point_vec);
    for (std::size_t k = 1; k <= number_of_subdivisions[0]; ++k)
    {
        for (std::size_t j = 1; j <= number_of_subdivisions[1]; ++j)
        {
            auto const offset = (j - 1) * (number_of_subdivisions[0] + 1);
            auto const upper_offset = j * (number_of_subdivisions[0] + 1);
            surface.addTriangle(offset + k - 1, offset + k,
                                upper_offset + k - 1);
            surface.addTriangle(offset + k, upper_offset + k,
                                upper_offset + k - 1);
        }
    }

    GeoLib::markUsedPoints(surface, used_points);

    ASSERT_EQ(point_vec.size(),
              std::count(used_points.begin(), used_points.end(), true));
}

TEST_F(GeoLibMarkUsedPointsInSurface, AllPointsExceptOneUsedInTriangles)
{
    // create surface covering all points except one corner point
    auto surface = GeoLib::Surface(point_vec);
    for (std::size_t k = 1; k <= number_of_subdivisions[0]; ++k)
    {
        for (std::size_t j = 1; j <= number_of_subdivisions[1]; ++j)
        {
            auto const offset = (j - 1) * (number_of_subdivisions[0] + 1);
            auto const upper_offset = j * (number_of_subdivisions[0] + 1);
            surface.addTriangle(offset + k - 1, offset + k,
                                upper_offset + k - 1);
        }
    }

    GeoLib::markUsedPoints(surface, used_points);

    ASSERT_EQ(point_vec.size() - 1,  // subtract corner point
              std::count(used_points.begin(), used_points.end(), true));
}
