/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"
#include <ctime>
#include <random>

#include "GeoLib/PointVec.h"

class PointVecTest : public testing::Test
{
public:
	typedef std::vector<GeoLib::Point*> VectorOfPoints;

	PointVecTest()
		: gen(std::random_device() ()), name("JustAName")
	{
		ps_ptr = new VectorOfPoints;
	}

protected:
	// Generates n new points according to given random number distribution,
	// which is uniform distribution in [-1, 1]^3.
	void
	generateRandomPoints(std::size_t const n = 1000)
	{
		std::uniform_real_distribution<double> rnd(-1, 1);
		std::generate_n(std::back_inserter(*ps_ptr), n,
			[&]() { return new GeoLib::Point(rnd(gen), rnd(gen), rnd(gen)); });
	}

protected:
	std::mt19937 gen;
	const std::string name;

	VectorOfPoints* ps_ptr;
};

// Testing nullptr input vector.
TEST_F(PointVecTest, TestPointVecCtorNullptr)
{
	ASSERT_ANY_THROW(GeoLib::PointVec(name, nullptr));
	delete ps_ptr;
}

// Testing empty input vector.
TEST_F(PointVecTest, TestPointVecCtorEmpty)
{
	ASSERT_ANY_THROW(GeoLib::PointVec(name, ps_ptr));
}

// Testing input vector with single point.
TEST_F(PointVecTest, TestPointVecCtorSinglePoint)
{
	ps_ptr->push_back(new GeoLib::Point(0,0,0));
	GeoLib::PointVec* point_vec = nullptr;
	ASSERT_NO_THROW(point_vec = new GeoLib::PointVec(name, ps_ptr));
	ASSERT_EQ(std::size_t(1), point_vec->size());

	delete point_vec;
}

// Testing input vector with two different points.
TEST_F(PointVecTest, TestPointVecCtorTwoDiffPoints)
{
	ps_ptr->push_back(new GeoLib::Point(0,0,0));
	ps_ptr->push_back(new GeoLib::Point(1,0,0));

	GeoLib::PointVec* point_vec = nullptr;
	ASSERT_NO_THROW(point_vec = new GeoLib::PointVec(name, ps_ptr));
	ASSERT_EQ(std::size_t(2), point_vec->size());

	delete point_vec;
}

// Testing input vector with two equal points.
TEST_F(PointVecTest, TestPointVecCtorTwoEqualPoints)
{
	ps_ptr->push_back(new GeoLib::Point(0,0,0));
	ps_ptr->push_back(new GeoLib::Point(0,0,0));

	GeoLib::PointVec* point_vec = nullptr;
	ASSERT_NO_THROW(point_vec = new GeoLib::PointVec(name, ps_ptr));
	ASSERT_EQ(std::size_t(1), point_vec->size());

	delete point_vec;
}

// Testing input vector with single point.
TEST_F(PointVecTest, TestPointVecPushBack)
{
	ps_ptr->push_back(new GeoLib::Point(0,0,0));
	ps_ptr->push_back(new GeoLib::Point(1,0,0));
	ps_ptr->push_back(new GeoLib::Point(0,1,0));
	ps_ptr->push_back(new GeoLib::Point(0,0,1));
	GeoLib::PointVec point_vec(name, ps_ptr);

	// Adding a new point with same coordinates changes nothing.
	point_vec.push_back(new GeoLib::Point(0,0,0));
	point_vec.push_back(new GeoLib::Point(1,0,0));
	point_vec.push_back(new GeoLib::Point(0,1,0));
	point_vec.push_back(new GeoLib::Point(0,0,1));

	ASSERT_EQ(std::size_t(4), point_vec.size());
}

// Testing random input points.
TEST_F(PointVecTest, TestPointVecCtorRandomPoints)
{
	generateRandomPoints(10000);

	GeoLib::PointVec* point_vec = nullptr;
	ASSERT_NO_THROW(point_vec = new GeoLib::PointVec(name, ps_ptr));

	delete point_vec;
}
TEST_F(PointVecTest, TestPointVecCtorRandomPointsLargeEps)
{
	generateRandomPoints(10000);

	GeoLib::PointVec* point_vec = nullptr;
	ASSERT_NO_THROW(point_vec = new GeoLib::PointVec(name, ps_ptr,
					nullptr, GeoLib::PointVec::PointType::POINT, 1e-2));

	delete point_vec;
}
