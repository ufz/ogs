/**
 * @file TestPoint.cpp
 * @author Thomas Fischer
 * @date Nov 8, 2012
 *
 * @copyright
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "Point.h"

using namespace GeoLib;

TEST(GeoLib, PointComparisonLessEq)
{
	// first coordinate
	ASSERT_TRUE(lessEq(Point(0,1,1),Point(1,1,1)));
	ASSERT_FALSE(lessEq(Point(1,1,1),Point(0,1,1)));
	// second coordinate
	ASSERT_TRUE(lessEq(Point(1,0,1),Point(1,1,1)));
	ASSERT_FALSE(lessEq(Point(1,1,1),Point(1,0,1)));
	// third coordinate
	ASSERT_TRUE(lessEq(Point(1,1,0),Point(1,1,1)));
	ASSERT_FALSE(lessEq(Point(1,1,1),Point(1,1,0)));

	const double e(2*std::numeric_limits<double>::epsilon());
	// first coordinate
	ASSERT_TRUE(lessEq(Point(1-e,1,1),Point(1,1,1)));
	ASSERT_FALSE(lessEq(Point(1,1,1),Point(1-e,1,1)));
	// second coordinate
	ASSERT_TRUE(lessEq(Point(1,1-e,1),Point(1,1,1)));
	ASSERT_FALSE(lessEq(Point(1,1,1),Point(1,1-e,1)));
	// third coordinate
	ASSERT_TRUE(lessEq(Point(1,1,1-e),Point(1,1,1)));
	ASSERT_FALSE(lessEq(Point(1,1,1),Point(1,1,1-e)));

	ASSERT_TRUE(lessEq(Point(1,1,1),Point(1,1,1)));

	const double half_eps(0.5*std::numeric_limits<double>::epsilon());
	ASSERT_TRUE(lessEq(Point(1+half_eps,1,1),Point(1,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1+half_eps,1),Point(1,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1+half_eps),Point(1,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1),Point(1+half_eps,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1),Point(1,1+half_eps,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1),Point(1,1,1+half_eps)));

	ASSERT_TRUE(lessEq(Point(1-half_eps,1,1),Point(1,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1-half_eps,1),Point(1,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1-half_eps),Point(1,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1),Point(1-half_eps,1,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1),Point(1,1-half_eps,1)));
	ASSERT_TRUE(lessEq(Point(1,1,1),Point(1,1,1-half_eps)));

	const double m(std::numeric_limits<double>::min());
	ASSERT_TRUE(lessEq(Point(m+half_eps,m,m),Point(m,m,m)));
	ASSERT_TRUE(lessEq(Point(m,m+half_eps,m),Point(m,m,m)));
	ASSERT_TRUE(lessEq(Point(m,m,m+half_eps),Point(m,m,m)));
	ASSERT_TRUE(lessEq(Point(m,m,m),Point(m+half_eps,m,m)));
	ASSERT_TRUE(lessEq(Point(m,m,m),Point(m,m+half_eps,m)));
	ASSERT_TRUE(lessEq(Point(m,m,m),Point(m,m,m+half_eps)));

	const double zero(0.0);
	ASSERT_TRUE(lessEq(Point(zero+half_eps,zero,zero),Point(zero,zero,zero)));
	ASSERT_TRUE(lessEq(Point(zero,zero+half_eps,zero),Point(zero,zero,zero)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero+half_eps),Point(zero,zero,zero)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero),Point(zero+half_eps,zero,zero)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero),Point(zero,zero+half_eps,zero)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero),Point(zero,zero,zero+half_eps)));

	ASSERT_TRUE(lessEq(Point(m+half_eps,m,m),Point(zero,zero,zero)));
	ASSERT_TRUE(lessEq(Point(m,m+half_eps,m),Point(zero,zero,zero)));
	ASSERT_TRUE(lessEq(Point(m,m,m+half_eps),Point(zero,zero,zero)));
	ASSERT_TRUE(lessEq(Point(m,m,m),Point(zero+half_eps,zero,zero)));
	ASSERT_TRUE(lessEq(Point(m,m,m),Point(zero,zero+half_eps,zero)));
	ASSERT_TRUE(lessEq(Point(m,m,m),Point(zero,zero,zero+half_eps)));

	ASSERT_TRUE(lessEq(Point(zero+half_eps,zero,zero),Point(m,m,m)));
	ASSERT_TRUE(lessEq(Point(zero,zero+half_eps,zero),Point(m,m,m)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero+half_eps),Point(m,m,m)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero),Point(m+half_eps,m,m)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero),Point(m,m+half_eps,m)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero),Point(m,m,m+half_eps)));

	ASSERT_TRUE(lessEq(Point(half_eps+half_eps,half_eps,half_eps),Point(half_eps,half_eps,half_eps)));
	ASSERT_TRUE(lessEq(Point(half_eps,half_eps+half_eps,zero),Point(half_eps,half_eps,half_eps)));
	ASSERT_TRUE(lessEq(Point(zero,zero,zero+half_eps),Point(half_eps,half_eps,half_eps)));
	ASSERT_TRUE(lessEq(Point(half_eps,half_eps,half_eps),Point(half_eps+half_eps,half_eps,half_eps)));
	ASSERT_TRUE(lessEq(Point(half_eps,half_eps,half_eps),Point(half_eps,half_eps+half_eps,half_eps)));
	ASSERT_TRUE(lessEq(Point(half_eps,half_eps,half_eps),Point(half_eps,half_eps,half_eps+half_eps)));

	ASSERT_TRUE(lessEq(Point(10.0+half_eps,10.0,10.0),Point(10.0,10.0,10.0)));
	ASSERT_TRUE(lessEq(Point(10.0,10.0+half_eps,10.0),Point(10.0,10.0,10.0)));
	ASSERT_TRUE(lessEq(Point(10.0,10.0,10.0+half_eps),Point(10.0,10.0,10.0)));
	ASSERT_TRUE(lessEq(Point(10.0,10.0,10.0),Point(10.0+half_eps,10.0,10.0)));
	ASSERT_TRUE(lessEq(Point(10.0,10.0,10.0),Point(10.0,10.0+half_eps,10.0)));
	ASSERT_TRUE(lessEq(Point(10.0,10.0,10.0),Point(10.0,10.0,10.0+half_eps)));
}

TEST(GeoLib, PointComparisonOperatorLessEq)
{
	const double my_eps(std::numeric_limits<double>::epsilon());
	ASSERT_FALSE(Point(1.0+my_eps,1.0,1.0) <= Point(1.0,1.0,1.0));
	ASSERT_FALSE(Point(1.0,1.0+my_eps,1.0) <= Point(1.0,1.0,1.0));
	ASSERT_FALSE(Point(1.0,1.0,1.0+my_eps) <= Point(1.0,1.0,1.0));
	ASSERT_TRUE(Point(1.0,1.0,1.0) <= Point(1.0+my_eps,1.0,1.0));
	ASSERT_TRUE(Point(1.0,1.0,1.0) <= Point(1.0,1.0+my_eps,1.0));
	ASSERT_TRUE(Point(1.0,1.0,1.0) <= Point(1.0,1.0,1.0+my_eps));

	ASSERT_TRUE(Point(1.0-my_eps,1.0,1.0) <= Point(1.0,1.0,1.0));
	ASSERT_TRUE(Point(1.0,1.0-my_eps,1.0) <= Point(1.0,1.0,1.0));
	ASSERT_TRUE(Point(1.0,1.0,1.0-my_eps) <= Point(1.0,1.0,1.0));
	ASSERT_FALSE(Point(1.0,1.0,1.0) <= Point(1.0-my_eps,1.0,1.0));
	ASSERT_FALSE(Point(1.0,1.0,1.0) <= Point(1.0,1.0-my_eps,1.0));
	ASSERT_FALSE(Point(1.0,1.0,1.0) <= Point(1.0,1.0,1.0-my_eps));

	std::size_t n(10000);
	srand ( time(NULL) );
	for (std::size_t k(0); k<n; ++k) {
		double random_val_x(((double)(rand()) / RAND_MAX - 0.5)); //real_dist(rng));
		double random_val_y(((double)(rand()) / RAND_MAX - 0.5)); //real_dist(rng));
		double random_val_z(((double)(rand()) / RAND_MAX - 0.5)); //real_dist(rng));

		double big_x(random_val_x * std::numeric_limits<double>::max());
		double big_y(random_val_y * std::numeric_limits<double>::max());
		double big_z(random_val_z * std::numeric_limits<double>::max());

		ASSERT_TRUE(Point(big_x-my_eps,big_y,big_z) <= Point(big_x,big_y,big_z));
		ASSERT_TRUE(Point(big_x,big_y-my_eps,big_z) <= Point(big_x,big_y,big_z));
		ASSERT_TRUE(Point(big_x,big_y,big_z-my_eps) <= Point(big_x,big_y,big_z));

		ASSERT_TRUE(Point(big_x-my_eps,big_y-my_eps,big_z) <= Point(big_x,big_y,big_z));
		ASSERT_TRUE(Point(big_x-my_eps,big_y,big_z-my_eps) <= Point(big_x,big_y,big_z));
		ASSERT_TRUE(Point(big_x,big_y-my_eps,big_z-my_eps) <= Point(big_x,big_y,big_z));

		ASSERT_TRUE(Point(big_x-my_eps,big_y-my_eps,big_z-my_eps) <= Point(big_x,big_y,big_z));

		double small_x(random_val_x * std::numeric_limits<double>::epsilon());
		double small_y(random_val_y * std::numeric_limits<double>::epsilon());
		double small_z(random_val_z * std::numeric_limits<double>::epsilon());

		ASSERT_TRUE(Point(small_x-my_eps,small_y,small_z) <= Point(small_x,small_y,small_z));
		ASSERT_TRUE(Point(small_x,small_y-my_eps,small_z) <= Point(small_x,small_y,small_z));
		ASSERT_TRUE(Point(small_x,small_y,small_z-my_eps) <= Point(small_x,small_y,small_z));

		ASSERT_TRUE(Point(small_x-my_eps,small_y-my_eps,small_z) <= Point(small_x,small_y,small_z));
		ASSERT_TRUE(Point(small_x-my_eps,small_y,small_z-my_eps) <= Point(small_x,small_y,small_z));
		ASSERT_TRUE(Point(small_x,small_y-my_eps,small_z-my_eps) <= Point(small_x,small_y,small_z));

		ASSERT_TRUE(Point(small_x-my_eps,small_y-my_eps,small_z-my_eps) <= Point(small_x,small_y,small_z));
	}
}
