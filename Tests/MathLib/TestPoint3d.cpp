/**
 * @file TestPoint3d.cpp
 * @author Thomas Fischer
 * @date Nov 8, 2012
 *
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "Point3d.h"
#include "MathTools.h"

using namespace MathLib;

TEST(MathLib, Point3dComparisonLessEq)
{
	// first coordinate
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{0,1,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_FALSE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{0,1,1}}))));
	// second coordinate
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,0,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_FALSE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,0,1}}))));
	// third coordinate
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,0}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_FALSE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1,0}}))));

	const double e(2*std::numeric_limits<double>::epsilon());
	// first coordinate
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1-e,1,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_FALSE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1-e,1,1}}))));
	// second coordinate
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1-e,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_FALSE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1-e,1}}))));
	// third coordinate
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1-e}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_FALSE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1,1-e}}))));

	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1,1}}))));

	const double half_eps(0.5*std::numeric_limits<double>::epsilon());
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1+half_eps,1,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1+half_eps,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1+half_eps}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1+half_eps,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1+half_eps,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1,1+half_eps}}))));

	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1-half_eps,1,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1-half_eps,1}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1-half_eps}})),Point3d(std::array<double,3>({{1,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1-half_eps,1,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1-half_eps,1}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{1,1,1}})),Point3d(std::array<double,3>({{1,1,1-half_eps}}))));

	const double m(std::numeric_limits<double>::min());
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m+half_eps,m,m}})),Point3d(std::array<double,3>({{m,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m+half_eps,m}})),Point3d(std::array<double,3>({{m,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m+half_eps}})),Point3d(std::array<double,3>({{m,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m}})),Point3d(std::array<double,3>({{m+half_eps,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m}})),Point3d(std::array<double,3>({{m,m+half_eps,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m}})),Point3d(std::array<double,3>({{m,m,m+half_eps}}))));

	const double zero(0.0);
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero+half_eps,zero,zero}})),Point3d(std::array<double,3>({{zero,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero+half_eps,zero}})),Point3d(std::array<double,3>({{zero,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero+half_eps}})),Point3d(std::array<double,3>({{zero,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero}})),Point3d(std::array<double,3>({{zero+half_eps,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero}})),Point3d(std::array<double,3>({{zero,zero+half_eps,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero}})),Point3d(std::array<double,3>({{zero,zero,zero+half_eps}}))));

	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m+half_eps,m,m}})),Point3d(std::array<double,3>({{zero,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m+half_eps,m}})),Point3d(std::array<double,3>({{zero,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m+half_eps}})),Point3d(std::array<double,3>({{zero,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m}})),Point3d(std::array<double,3>({{zero+half_eps,zero,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m}})),Point3d(std::array<double,3>({{zero,zero+half_eps,zero}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{m,m,m}})),Point3d(std::array<double,3>({{zero,zero,zero+half_eps}}))));

	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero+half_eps,zero,zero}})),Point3d(std::array<double,3>({{m,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero+half_eps,zero}})),Point3d(std::array<double,3>({{m,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero+half_eps}})),Point3d(std::array<double,3>({{m,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero}})),Point3d(std::array<double,3>({{m+half_eps,m,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero}})),Point3d(std::array<double,3>({{m,m+half_eps,m}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero}})),Point3d(std::array<double,3>({{m,m,m+half_eps}}))));

	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{half_eps+half_eps,half_eps,half_eps}})),Point3d(std::array<double,3>({{half_eps,half_eps,half_eps}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{half_eps,half_eps+half_eps,zero}})),Point3d(std::array<double,3>({{half_eps,half_eps,half_eps}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{zero,zero,zero+half_eps}})),Point3d(std::array<double,3>({{half_eps,half_eps,half_eps}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{half_eps,half_eps,half_eps}})),Point3d(std::array<double,3>({{half_eps+half_eps,half_eps,half_eps}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{half_eps,half_eps,half_eps}})),Point3d(std::array<double,3>({{half_eps,half_eps+half_eps,half_eps}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{half_eps,half_eps,half_eps}})),Point3d(std::array<double,3>({{half_eps,half_eps,half_eps+half_eps}}))));

	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{10.0+half_eps,10.0,10.0}})),Point3d(std::array<double,3>({{10.0,10.0,10.0}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{10.0,10.0+half_eps,10.0}})),Point3d(std::array<double,3>({{10.0,10.0,10.0}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{10.0,10.0,10.0+half_eps}})),Point3d(std::array<double,3>({{10.0,10.0,10.0}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{10.0,10.0,10.0}})),Point3d(std::array<double,3>({{10.0+half_eps,10.0,10.0}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{10.0,10.0,10.0}})),Point3d(std::array<double,3>({{10.0,10.0+half_eps,10.0}}))));
	ASSERT_TRUE(lessEq(Point3d(std::array<double,3>({{10.0,10.0,10.0}})),Point3d(std::array<double,3>({{10.0,10.0,10.0+half_eps}}))));
}

TEST(MathLib, Point3dComparisonOperatorLessEq)
{
	const double my_eps(std::numeric_limits<double>::epsilon());
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0+my_eps,1.0,1.0}})) <=
	Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0+my_eps,1.0}})) <=
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,1.0+my_eps}})) <=
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <=
		Point3d(std::array<double,3>({{1.0+my_eps,1.0,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <=
		Point3d(std::array<double,3>({{1.0,1.0+my_eps,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <=
		Point3d(std::array<double,3>({{1.0,1.0,1.0+my_eps}})));

	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0-my_eps,1.0,1.0}})) <=
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0-my_eps,1.0}})) <=
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,1.0-my_eps}})) <=
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <=
		Point3d(std::array<double,3>({{1.0-my_eps,1.0,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <=
		Point3d(std::array<double,3>({{1.0,1.0-my_eps,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <=
		Point3d(std::array<double,3>({{1.0,1.0,1.0-my_eps}})));

	std::size_t n(10000);
	srand ( static_cast<unsigned>(time(NULL)) );
	for (std::size_t k(0); k<n; ++k) {
		double random_val_x(((double)(rand()) / RAND_MAX - 0.5)); //real_dist(rng));
		double random_val_y(((double)(rand()) / RAND_MAX - 0.5)); //real_dist(rng));
		double random_val_z(((double)(rand()) / RAND_MAX - 0.5)); //real_dist(rng));

		double big_x(random_val_x * std::numeric_limits<double>::max());
		double big_y(random_val_y * std::numeric_limits<double>::max());
		double big_z(random_val_z * std::numeric_limits<double>::max());

		ASSERT_TRUE(Point3d(std::array<double,3>({{big_x-my_eps,big_y,big_z}})) <= Point3d(std::array<double,3>({{big_x,big_y,big_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{big_x,big_y-my_eps,big_z}})) <= Point3d(std::array<double,3>({{big_x,big_y,big_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{big_x,big_y,big_z-my_eps}})) <= Point3d(std::array<double,3>({{big_x,big_y,big_z}})));

		ASSERT_TRUE(Point3d(std::array<double,3>({{big_x-my_eps,big_y-my_eps,big_z}})) <= Point3d(std::array<double,3>({{big_x,big_y,big_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{big_x-my_eps,big_y,big_z-my_eps}})) <= Point3d(std::array<double,3>({{big_x,big_y,big_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{big_x,big_y-my_eps,big_z-my_eps}})) <= Point3d(std::array<double,3>({{big_x,big_y,big_z}})));

		ASSERT_TRUE(Point3d(std::array<double,3>({{big_x-my_eps,big_y-my_eps,big_z-my_eps}})) <= Point3d(std::array<double,3>({{big_x,big_y,big_z}})));

		double small_x(random_val_x * std::numeric_limits<double>::epsilon());
		double small_y(random_val_y * std::numeric_limits<double>::epsilon());
		double small_z(random_val_z * std::numeric_limits<double>::epsilon());

		ASSERT_TRUE(Point3d(std::array<double,3>({{small_x-my_eps,small_y,small_z}})) <= Point3d(std::array<double,3>({{small_x,small_y,small_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{small_x,small_y-my_eps,small_z}})) <= Point3d(std::array<double,3>({{small_x,small_y,small_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{small_x,small_y,small_z-my_eps}})) <= Point3d(std::array<double,3>({{small_x,small_y,small_z}})));

		ASSERT_TRUE(Point3d(std::array<double,3>({{small_x-my_eps,small_y-my_eps,small_z}})) <= Point3d(std::array<double,3>({{small_x,small_y,small_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{small_x-my_eps,small_y,small_z-my_eps}})) <= Point3d(std::array<double,3>({{small_x,small_y,small_z}})));
		ASSERT_TRUE(Point3d(std::array<double,3>({{small_x,small_y-my_eps,small_z-my_eps}})) <= Point3d(std::array<double,3>({{small_x,small_y,small_z}})));

		ASSERT_TRUE(Point3d(std::array<double,3>({{small_x-my_eps,small_y-my_eps,small_z-my_eps}})) <= Point3d(std::array<double,3>({{small_x,small_y,small_z}})));
	}
}

// test for operator==
TEST(MathLib, Point3dComparisonOperatorEqual)
{
	srand(static_cast<unsigned>(time(nullptr)));
	double x0(((double)(rand()) / RAND_MAX - 0.5));
	double x1(((double)(rand()) / RAND_MAX - 0.5));
	double x2(((double)(rand()) / RAND_MAX - 0.5));

	MathLib::Point3d a(std::array<double,3>({{x0, x1, x2}}));
	MathLib::Point3d b(std::array<double,3>({{x0, x1, x2}}));
	ASSERT_TRUE(((a <= b) && (b <= a)) == (a == b));
	ASSERT_TRUE(a == b);
	ASSERT_TRUE((lessEq(a,b) && lessEq(b,a)) == (a == b));

	double tol(std::numeric_limits<double>::min());
	b[2] += tol;
	ASSERT_TRUE(((a <= b) && (b <= a)) == (a == b));
	b[1] = 0.0;
	b[2] = 0.0;
	ASSERT_TRUE(((a <= b) && (b <= a)) == (a == b));

	tol = std::numeric_limits<double>::epsilon();
	ASSERT_FALSE(Point3d(std::array<double,3>({{tol,1.0,1.0}})) == Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,tol,1.0}})) == Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,tol}})) == Point3d(std::array<double,3>({{1.0,1.0,1.0}})));

	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) == Point3d(std::array<double,3>({{1.0+tol,1.0,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) == Point3d(std::array<double,3>({{1.0,1.0+tol,1.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) == Point3d(std::array<double,3>({{1.0,1.0,1.0+tol}})));

	// very small difference in one coordinate
	tol = std::numeric_limits<double>::min();
	ASSERT_TRUE(Point3d(std::array<double,3>({{tol,0.0,0.0}})) == Point3d(std::array<double,3>({{0.0,0.0,0.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,tol,0.0}})) == Point3d(std::array<double,3>({{0.0,0.0,0.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,0.0,tol}})) == Point3d(std::array<double,3>({{0.0,0.0,0.0}})));

	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,0.0,0.0}})) == Point3d(std::array<double,3>({{tol,0.0,0.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,0.0,0.0}})) == Point3d(std::array<double,3>({{0.0,tol,0.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,0.0,0.0}})) == Point3d(std::array<double,3>({{0.0,0.0,tol}})));

	a = Point3d(std::array<double,3>({{0.0,0.0,0.0}}));
	b = Point3d(std::array<double,3>({{0.0,0.0,0.0}}));
	a[0] = pow(std::numeric_limits<double>::epsilon(),5);
	ASSERT_TRUE((lessEq(a,b) && lessEq(b,a)) == (a == b));
}


// test for operator<
TEST(MathLib, Point3dComparisonOperatorLess)
{
	srand(static_cast<unsigned>(time(nullptr)));
	double x0(((double)(rand()) / RAND_MAX - 0.5));
	double x1(((double)(rand()) / RAND_MAX - 0.5));
	double x2(((double)(rand()) / RAND_MAX - 0.5));

	MathLib::Point3d a(std::array<double,3>({{x0, x1, x2}}));
	MathLib::Point3d b(std::array<double,3>({{x0, x1, x2}}));
	ASSERT_FALSE((a < b) && (b < a));

	double tol = std::numeric_limits<double>::epsilon();
	ASSERT_TRUE(Point3d(std::array<double,3>({{tol,1.0,1.0}})) <
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,tol,1.0}})) <
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,tol}})) <
		Point3d(std::array<double,3>({{1.0,1.0,1.0}})));

	// very small difference in one coordinate
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <
		Point3d(std::array<double,3>({{1.0+tol,1.0,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <
		Point3d(std::array<double,3>({{1.0,1.0+tol,1.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{1.0,1.0,1.0}})) <
		Point3d(std::array<double,3>({{1.0,1.0,1.0+tol}})));

	tol = std::numeric_limits<double>::min();
	ASSERT_FALSE(Point3d(std::array<double,3>({{tol,0.0,0.0}})) <
		Point3d(std::array<double,3>({{0.0,0.0,0.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{0.0,tol,0.0}})) <
		Point3d(std::array<double,3>({{0.0,0.0,0.0}})));
	ASSERT_FALSE(Point3d(std::array<double,3>({{0.0,0.0,tol}})) <
		Point3d(std::array<double,3>({{0.0,0.0,0.0}})));

	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,0.0,0.0}})) <
		Point3d(std::array<double,3>({{tol,0.0,0.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,0.0,0.0}})) <
		Point3d(std::array<double,3>({{0.0,tol,0.0}})));
	ASSERT_TRUE(Point3d(std::array<double,3>({{0.0,0.0,0.0}})) <
		Point3d(std::array<double,3>({{0.0,0.0,tol}})));
}

