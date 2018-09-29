/**
 * @file TestVector3.cpp
 * @date Feb 28, 2014
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <array>

#include "gtest/gtest.h"

#include "MathLib/Vector3.h"
#include "GeoLib/Point.h"

using namespace MathLib;

TEST(MathLib, TestVector3Constructor)
{
    // *** test default constructor
    Vector3 u;
    // test coordinates of default constructed vec
    ASSERT_NEAR(0.0, u[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.0, u[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.0, u[2], std::numeric_limits<double>::epsilon());

    // *** test constructor taking 3 double values
    Vector3 v(1.0, 3.0, 5.0);
    ASSERT_NEAR(1.0, v[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(3.0, v[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(5.0, v[2], std::numeric_limits<double>::epsilon());

    // *** test copy constructor
    Vector3 v_copy(v);
    // test equality of coordinates
    ASSERT_NEAR(v[0], v_copy[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(v[1], v_copy[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(v[2], v_copy[2], std::numeric_limits<double>::epsilon());

    // *** test constructor taking TemplatePoint
    std::array<double,3> ap = {{0, 1, 2}};
    TemplatePoint<double> p(ap);
    Vector3 vp(p);
    ASSERT_NEAR(0.0, vp[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(1.0, vp[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(2.0, vp[2], std::numeric_limits<double>::epsilon());

    // *** test constructing Vector from two TemplatePoints
    std::array<double,3> aa = {{1, 2, 3}}; // necessary for old compilers
    std::array<double,3> ab = {{6, 5, 4}}; // necessary for old compilers
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
    ASSERT_NEAR(1.0, v[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(3.0, v[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(5.0, v[2], std::numeric_limits<double>::epsilon());

    Vector3 w(5.0, 3.0, 1.0);
    // operator+
    Vector3 res(v+w);
    ASSERT_NEAR(6.0, res[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(6.0, res[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(6.0, res[2], std::numeric_limits<double>::epsilon());

    // operator-
    res = v-w;
    ASSERT_NEAR(-4.0, res[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR( 0.0, res[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR( 4.0, res[2], std::numeric_limits<double>::epsilon());

    // test operator*=
    v *= 2.0;
    ASSERT_NEAR(2.0, v[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(6.0, v[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(10.0, v[2], std::numeric_limits<double>::epsilon());

    // test operator+=
    v += w;
    ASSERT_NEAR(7.0, v[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(9.0, v[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(11.0, v[2], std::numeric_limits<double>::epsilon());

    // test operator-=
    v -= w;
    ASSERT_NEAR(2.0, v[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(6.0, v[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(10.0, v[2], std::numeric_limits<double>::epsilon());
}

TEST(MathLib, TestVector3Multiplications)
{
    // test scalar product
    Vector3 v(1.0, 3.0, 5.0);
    Vector3 w(3.0, -2.0, 1.0);

    ASSERT_NEAR(2.0, scalarProduct(v,w), std::numeric_limits<double>::epsilon());

    // test cross product
    Vector3 e1(1.0, 0.0, 0.0);
    Vector3 e2(0.0, 1.0, 0.0);
    Vector3 e3(0.0, 0.0, 1.0);

    Vector3 res_e1e2(crossProduct(e1, e2)); // should be e3
    ASSERT_NEAR(e3[0], res_e1e2[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(e3[1], res_e1e2[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(e3[2], res_e1e2[2], std::numeric_limits<double>::epsilon());

    Vector3 res_e2e3(crossProduct(e2, e3)); // should be e1
    ASSERT_NEAR(e1[0], res_e2e3[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(e1[1], res_e2e3[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(e1[2], res_e2e3[2], std::numeric_limits<double>::epsilon());

    Vector3 res_e3e1(crossProduct(e3, e1)); // should be e2
    ASSERT_NEAR(e2[0], res_e3e1[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(e2[1], res_e3e1[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(e2[2], res_e3e1[2], std::numeric_limits<double>::epsilon());

    Vector3 res_e2e1(crossProduct(e2, e1)); // should be -e3
    ASSERT_NEAR(-e3[0], res_e2e1[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-e3[1], res_e2e1[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-e3[2], res_e2e1[2], std::numeric_limits<double>::epsilon());

    Vector3 res_e3e2(crossProduct(e3, e2)); // should be -e1
    ASSERT_NEAR(-e1[0], res_e3e2[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-e1[1], res_e3e2[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-e1[2], res_e3e2[2], std::numeric_limits<double>::epsilon());

    Vector3 res_e1e3(crossProduct(e1, e3)); // should be -e2
    ASSERT_NEAR(-e2[0], res_e1e3[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-e2[1], res_e1e3[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-e2[2], res_e1e3[2], std::numeric_limits<double>::epsilon());

    // test multplication with scalar
    v = -1.0 * v;
    ASSERT_NEAR(-1.0, v[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-3.0, v[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-5.0, v[2], std::numeric_limits<double>::epsilon());
    v = v * -1.0;
    ASSERT_NEAR(1.0, v[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(3.0, v[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(5.0, v[2], std::numeric_limits<double>::epsilon());

    // test normalisation
    v.normalize();
    ASSERT_NEAR(1.0, v.getLength(), std::numeric_limits<double>::epsilon());

}

