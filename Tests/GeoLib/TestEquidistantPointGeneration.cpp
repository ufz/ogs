/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <boost/math/constants/constants.hpp>
#include <random>

#include "GeoLib/Utils.h"

struct GeoLibGenerateEquidistantPoints : public testing::Test
{
    GeoLibGenerateEquidistantPoints()
    {
        // enable random engine
        std::random_device rd;
        std::mt19937 random_engine_mt19937(rd());
        constexpr auto pi = boost::math::double_constants::pi;
        std::normal_distribution<> normal_dist_phi(-pi, pi);  // azimuthal angle
        std::normal_distribution<> normal_dist_theta(0, pi);  // polar angle
        // generate point on unit sphere
        double const phi = normal_dist_phi(random_engine_mt19937);
        double const theta = normal_dist_theta(random_engine_mt19937);

        start = MathLib::Point3d{{std::cos(phi) * std::sin(theta),
                                  std::sin(phi) * std::sin(theta),
                                  std::cos(theta)}};
        // mirror of start point at origin
        end = MathLib::Point3d{{-start[0], -start[1], -start[2]}};
    }

    ~GeoLibGenerateEquidistantPoints()
    {
        for (auto* p : equidistant_points)
        {
            delete p;
        }
    }

    void checkLength() const
    {
        constexpr double eps = 4 * std::numeric_limits<double>::epsilon();
        double sum = 0.0;
        for (std::size_t i = 0; i < equidistant_points.size() - 1; ++i)
        {
            sum += std::sqrt(MathLib::sqrDist(*(equidistant_points[i]),
                                              *(equidistant_points[i + 1])));
        }
        auto const sum_start_end = std::sqrt(MathLib::sqrDist(start, end));
        EXPECT_NEAR(sum_start_end, sum, eps)
            << start[0] << " " << start[1] << " " << start[2] << " -- "
            << end[0] << " " << end[1] << " " << end[2];
    }

    MathLib::Point3d start;
    MathLib::Point3d end;
    std::vector<GeoLib::Point*> equidistant_points;
};

TEST_F(GeoLibGenerateEquidistantPoints, IdenticalStartEndPoints)
{
    equidistant_points = GeoLib::generateEquidistantPoints(start, start, 2);
    ASSERT_EQ(4, equidistant_points.size());
    EXPECT_DOUBLE_EQ(start[0], (*equidistant_points[0])[0]);
    EXPECT_DOUBLE_EQ(start[1], (*equidistant_points[0])[1]);
    EXPECT_DOUBLE_EQ(start[2], (*equidistant_points[0])[2]);
    EXPECT_DOUBLE_EQ(start[0], (*equidistant_points[1])[0]);
    EXPECT_DOUBLE_EQ(start[1], (*equidistant_points[1])[1]);
    EXPECT_DOUBLE_EQ(start[2], (*equidistant_points[1])[2]);
}

TEST_F(GeoLibGenerateEquidistantPoints, ZeroSubdivisions)
{
    equidistant_points = GeoLib::generateEquidistantPoints(start, end, 0);
    ASSERT_EQ(2, equidistant_points.size());
    checkLength();
}

TEST_F(GeoLibGenerateEquidistantPoints, OneSubdivision)
{
    equidistant_points = GeoLib::generateEquidistantPoints(start, end, 1);
    ASSERT_EQ(3, equidistant_points.size());
    EXPECT_DOUBLE_EQ(0.0, (*equidistant_points[1])[0]);
    EXPECT_DOUBLE_EQ(0.0, (*equidistant_points[1])[1]);
    EXPECT_DOUBLE_EQ(0.0, (*equidistant_points[1])[2]);
    checkLength();
}

TEST_F(GeoLibGenerateEquidistantPoints, ThreeSubdivisions)
{
    equidistant_points = GeoLib::generateEquidistantPoints(start, end, 3);
    ASSERT_EQ(5, equidistant_points.size());
    EXPECT_DOUBLE_EQ(0.0, (*equidistant_points[2])[0]);
    EXPECT_DOUBLE_EQ(0.0, (*equidistant_points[2])[1]);
    EXPECT_DOUBLE_EQ(0.0, (*equidistant_points[2])[2]);
    checkLength();
}

TEST_F(GeoLibGenerateEquidistantPoints, FiveSubdivisions)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();
    equidistant_points = GeoLib::generateEquidistantPoints(start, end, 5);
    ASSERT_EQ(7, equidistant_points.size());
    EXPECT_NEAR(0.0, (*equidistant_points[3])[0], eps);
    EXPECT_NEAR(0.0, (*equidistant_points[3])[1], eps);
    EXPECT_NEAR(0.0, (*equidistant_points[3])[2], eps);
    EXPECT_NEAR(
        MathLib::sqrDist(start, *(equidistant_points[1])),
        MathLib::sqrDist(*(equidistant_points[1]), *(equidistant_points[2])),
        eps);
    EXPECT_NEAR(MathLib::sqrDist(start, *(equidistant_points[1])),
                MathLib::sqrDist(*(equidistant_points[5]), end),
                eps);
    checkLength();
}
