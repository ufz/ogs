/**
 * @file TestTemplatePoint.cpp
 * @date Mar 13, 2014
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <array>

#include "gtest/gtest.h"

#include "MathLib/TemplatePoint.h"

using namespace MathLib;

TEST(MathLib, TemplatePointConstructors)
{
    // *** test default constructors for different DIM values
    TemplatePoint<double,3> p0;
    // test coordinates of default constructed point
    ASSERT_EQ(0.0, p0[0]);
    ASSERT_EQ(0.0, p0[1]);
    ASSERT_EQ(0.0, p0[2]);

    TemplatePoint<double,2> p1;
    ASSERT_EQ(0.0, p1[0]);
    ASSERT_EQ(0.0, p1[1]);

    TemplatePoint<double,5> p2;
    ASSERT_EQ(0.0, p2[0]);
    ASSERT_EQ(0.0, p2[1]);
    ASSERT_EQ(0.0, p2[2]);
    ASSERT_EQ(0.0, p2[3]);
    ASSERT_EQ(0.0, p2[4]);

    // *** test copy constructor
    TemplatePoint<double,3> p0_copy(p0);
    // test equality of coordinates
    ASSERT_EQ(p0[0], p0_copy[0]);
    ASSERT_EQ(p0[1], p0_copy[1]);
    ASSERT_EQ(p0[2], p0_copy[2]);

    // *** test constructor taking std::array
    std::array<double,3> array = {{0, 1, 2}};
    TemplatePoint<double> p3(array);
    ASSERT_EQ(0.0, p3[0]);
    ASSERT_EQ(1.0, p3[1]);
    ASSERT_EQ(2.0, p3[2]);
}

TEST(MathLib, TemplatePointOperators)
{
    TemplatePoint<double,3> p;
    // access operator
    p[0] = 1.0;
    p[1] = 3.0;
    p[2] = 5.0;
    ASSERT_EQ(1.0, p[0]);
    ASSERT_EQ(3.0, p[1]);
    ASSERT_EQ(5.0, p[2]);
}

