/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <ctime>
#include <random>

#include "GeoLib/PointVec.h"

class PointVecTest : public testing::Test
{
public:
    using VectorOfPoints = std::vector<GeoLib::Point*>;

    PointVecTest() : gen(std::random_device()()), name("JustAName") {}

protected:
    // Generates n new points according to given random number distribution,
    // which is uniform distribution in [-1, 1]^3.
    void generateRandomPoints(VectorOfPoints& ps, std::size_t const n = 1000)
    {
        std::uniform_real_distribution<double> rnd(-1, 1);
        std::generate_n(std::back_inserter(ps), n,
                        [&]() {
                            return new GeoLib::Point(rnd(gen), rnd(gen),
                                                     rnd(gen), ps.size());
                        });
    }

protected:
    std::mt19937 gen;
    std::string const name;
};

// Testing empty input vector.
TEST_F(PointVecTest, TestPointVecCtorEmpty)
{
    VectorOfPoints points;
    ASSERT_THROW(GeoLib::PointVec(name, std::move(points)), std::runtime_error);
}

// Testing input vector with single point.
TEST_F(PointVecTest, TestPointVecCtorSinglePoint)
{
    VectorOfPoints points;
    points.push_back(new GeoLib::Point(0, 0, 0, 0));
    GeoLib::PointVec const point_vec(name, std::move(points));
    ASSERT_EQ(std::size_t(1), point_vec.size());
}

// Testing input vector with two different points.
TEST_F(PointVecTest, TestPointVecCtorTwoDiffPoints)
{
    VectorOfPoints points;
    points.push_back(new GeoLib::Point(0, 0, 0, 0));
    points.push_back(new GeoLib::Point(1, 0, 0, 1));

    GeoLib::PointVec const point_vec(name, std::move(points));
    ASSERT_EQ(std::size_t(2), point_vec.size());
}

// Testing input vector with two equal points.
TEST_F(PointVecTest, TestPointVecCtorTwoEqualPoints)
{
    VectorOfPoints points;
    points.push_back(new GeoLib::Point(0, 0, 0, 0));
    points.push_back(new GeoLib::Point(0, 0, 0, 1));

    GeoLib::PointVec const point_vec(name, std::move(points));
    ASSERT_EQ(std::size_t(1), point_vec.size());
}

// Testing input vector with single point.
TEST_F(PointVecTest, TestPointVecPushBack)
{
    VectorOfPoints points;
    points.push_back(new GeoLib::Point(0, 0, 0, 0));
    points.push_back(new GeoLib::Point(1, 0, 0, 1));
    points.push_back(new GeoLib::Point(0, 1, 0, 2));
    points.push_back(new GeoLib::Point(0, 0, 1, 3));

    GeoLib::PointVec point_vec(name, std::move(points));
    ASSERT_EQ(std::size_t(0), point_vec.getIDMap()[0]);
    ASSERT_EQ(std::size_t(1), point_vec.getIDMap()[1]);
    ASSERT_EQ(std::size_t(2), point_vec.getIDMap()[2]);
    ASSERT_EQ(std::size_t(3), point_vec.getIDMap()[3]);

    // Adding some points that are already existing.
    ASSERT_EQ(std::size_t(0),
              point_vec.push_back(new GeoLib::Point(0, 0, 0, 4)));
    ASSERT_EQ(std::size_t(1),
              point_vec.push_back(new GeoLib::Point(1, 0, 0, 5)));
    ASSERT_EQ(std::size_t(2),
              point_vec.push_back(new GeoLib::Point(0, 1, 0, 6)));
    ASSERT_EQ(std::size_t(3),
              point_vec.push_back(new GeoLib::Point(0, 0, 1, 7)));

    ASSERT_EQ(std::size_t(4), point_vec.size());
    ASSERT_EQ(std::size_t(8), point_vec.getIDMap().size());

    ASSERT_EQ(std::size_t(0), point_vec.getIDMap()[4]);
    ASSERT_EQ(std::size_t(1), point_vec.getIDMap()[5]);
    ASSERT_EQ(std::size_t(2), point_vec.getIDMap()[6]);
    ASSERT_EQ(std::size_t(3), point_vec.getIDMap()[7]);

    // Adding again some already existing points.
    ASSERT_EQ(std::size_t(0),
              point_vec.push_back(new GeoLib::Point(0, 0, 0, 8)));
    ASSERT_EQ(std::size_t(1),
              point_vec.push_back(new GeoLib::Point(1, 0, 0, 9)));
    ASSERT_EQ(std::size_t(2),
              point_vec.push_back(new GeoLib::Point(0, 1, 0, 10)));
    ASSERT_EQ(std::size_t(3),
              point_vec.push_back(new GeoLib::Point(0, 0, 1, 11)));

    ASSERT_EQ(std::size_t(4), point_vec.size());
    ASSERT_EQ(std::size_t(12), point_vec.getIDMap().size());

    ASSERT_EQ(std::size_t(0), point_vec.getIDMap()[8]);
    ASSERT_EQ(std::size_t(1), point_vec.getIDMap()[9]);
    ASSERT_EQ(std::size_t(2), point_vec.getIDMap()[10]);
    ASSERT_EQ(std::size_t(3), point_vec.getIDMap()[11]);

    // Adding some new points.
    ASSERT_EQ(std::size_t(4),
              point_vec.push_back(new GeoLib::Point(0.1, 0.1, 0.1, 12)));
    ASSERT_EQ(std::size_t(5),
              point_vec.push_back(new GeoLib::Point(1.1, 0.1, 0.1, 13)));
    ASSERT_EQ(std::size_t(6),
              point_vec.push_back(new GeoLib::Point(0.1, 1.1, 0.1, 14)));
    ASSERT_EQ(std::size_t(7),
              point_vec.push_back(new GeoLib::Point(0.1, 0.1, 1.1, 15)));

    ASSERT_EQ(std::size_t(8), point_vec.size());
    ASSERT_EQ(std::size_t(16), point_vec.getIDMap().size());

    ASSERT_EQ(std::size_t(4), point_vec.getIDMap()[12]);
    ASSERT_EQ(std::size_t(5), point_vec.getIDMap()[13]);
    ASSERT_EQ(std::size_t(6), point_vec.getIDMap()[14]);
    ASSERT_EQ(std::size_t(7), point_vec.getIDMap()[15]);

    // Adding again some already existing points.
    ASSERT_EQ(std::size_t(0),
              point_vec.push_back(new GeoLib::Point(0, 0, 0, 16)));
    ASSERT_EQ(std::size_t(1),
              point_vec.push_back(new GeoLib::Point(1, 0, 0, 17)));
    ASSERT_EQ(std::size_t(2),
              point_vec.push_back(new GeoLib::Point(0, 1, 0, 18)));
    ASSERT_EQ(std::size_t(3),
              point_vec.push_back(new GeoLib::Point(0, 0, 1, 19)));

    ASSERT_EQ(std::size_t(8), point_vec.size());
    ASSERT_EQ(std::size_t(20), point_vec.getIDMap().size());

    ASSERT_EQ(std::size_t(0), point_vec.getIDMap()[16]);
    ASSERT_EQ(std::size_t(1), point_vec.getIDMap()[17]);
    ASSERT_EQ(std::size_t(2), point_vec.getIDMap()[18]);
    ASSERT_EQ(std::size_t(3), point_vec.getIDMap()[19]);

    ASSERT_EQ(std::size_t(0), point_vec.getIDMap()[0]);
    ASSERT_EQ(std::size_t(1), point_vec.getIDMap()[1]);
    ASSERT_EQ(std::size_t(2), point_vec.getIDMap()[2]);
    ASSERT_EQ(std::size_t(3), point_vec.getIDMap()[3]);
    ASSERT_EQ(std::size_t(0), point_vec.getIDMap()[4]);
    ASSERT_EQ(std::size_t(1), point_vec.getIDMap()[5]);
    ASSERT_EQ(std::size_t(2), point_vec.getIDMap()[6]);
    ASSERT_EQ(std::size_t(3), point_vec.getIDMap()[7]);
    ASSERT_EQ(std::size_t(0), point_vec.getIDMap()[8]);
    ASSERT_EQ(std::size_t(1), point_vec.getIDMap()[9]);
    ASSERT_EQ(std::size_t(2), point_vec.getIDMap()[10]);
    ASSERT_EQ(std::size_t(3), point_vec.getIDMap()[11]);
    ASSERT_EQ(std::size_t(4), point_vec.getIDMap()[12]);
    ASSERT_EQ(std::size_t(5), point_vec.getIDMap()[13]);
    ASSERT_EQ(std::size_t(6), point_vec.getIDMap()[14]);
    ASSERT_EQ(std::size_t(7), point_vec.getIDMap()[15]);
}

// Testing random input points.
TEST_F(PointVecTest, TestPointVecCtorRandomPoints)
{
    VectorOfPoints points;
    generateRandomPoints(points, 10000);

    auto point_vec =
        std::make_unique<GeoLib::PointVec>(name, std::move(points));
}

TEST_F(PointVecTest, TestPointVecCtorRandomPointsLargeEps)
{
    VectorOfPoints points;
    generateRandomPoints(points, 10000);

    auto point_vec = std::make_unique<GeoLib::PointVec>(
        name, std::move(points), GeoLib::PointVec::PointType::POINT, 1e-2);
}
