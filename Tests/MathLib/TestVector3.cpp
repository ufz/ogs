/**
 * @file TestVector3.cpp
 * @date Feb 28, 2014
 *
 * @copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <array>

#include "gtest/gtest.h"

#include "Point.h"
#include "Vector3.h"

using namespace MathLib;

TEST(MathLib, TestVector3Constructor)
{
	// *** test default constructor
	Vector3 u;
	// test coordinates of default constructed vec
	ASSERT_NEAR(0.0, u[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(0.0, u[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(0.0, u[2], std::numeric_limits<double>::min());

	// *** test constructor taking 3 double values
	Vector3 v(1.0, 3.0, 5.0);
	ASSERT_NEAR(1.0, v[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(3.0, v[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(5.0, v[2], std::numeric_limits<double>::min());

	// *** test copy constructor
	Vector3 v_copy(v);
	// test equality of coordinates
	ASSERT_NEAR(v[0], v_copy[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(v[1], v_copy[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(v[2], v_copy[2], std::numeric_limits<double>::min());

	// *** test constructor taking TemplatePoint
	std::array<double,3> ap = {0, 1, 2};
	TemplatePoint<double> p(ap);
	Vector3 vp(p);
	ASSERT_NEAR(0.0, vp[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(1.0, vp[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(2.0, vp[2], std::numeric_limits<double>::min());

	// *** test constructing Vector from two TemplatePoints
	std::array<double,3> aa = {1, 2, 3}; // necessary for old compilers
	std::array<double,3> ab = {6, 5, 4}; // necessary for old compilers
	TemplatePoint<double,3> a(aa);
	TemplatePoint<double,3> b(ab);
	Vector3 w(a,b);
	// test coordinates of constructed Vector3 w = (b-a)
	ASSERT_NEAR(5.0, w[0], std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(3.0, w[1], std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(1.0, w[2], std::numeric_limits<double>::epsilon());
}

TEST(MathLib, TestVector3Operators)
{
	Vector3 v;
	// access operator
	v[0] = 1.0;
	v[1] = 3.0;
	v[2] = 5.0;
	ASSERT_NEAR(1.0, v[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(3.0, v[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(5.0, v[2], std::numeric_limits<double>::min());

	Vector3 w(5.0, 3.0, 1.0);
	// operator+
	Vector3 res(v+w);
	ASSERT_NEAR(6.0, res[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(6.0, res[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(6.0, res[2], std::numeric_limits<double>::min());

	// operator-
	res = v-w;
	ASSERT_NEAR(-4.0, res[0], std::numeric_limits<double>::min());
	ASSERT_NEAR( 0.0, res[1], std::numeric_limits<double>::min());
	ASSERT_NEAR( 4.0, res[2], std::numeric_limits<double>::min());
}
