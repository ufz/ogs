/**
 * @file TestWeightedPoint.cpp
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "MathLib/TemplateWeightedPoint.h"

TEST(MathLib, WeightedPoint1D)
{
    std::array<double, 1> pnt;
    pnt[0] = 0.5;
    double w = 100.0;
    MathLib::WeightedPoint1D wpnt_1d(pnt, w);

    ASSERT_NEAR(pnt[0], wpnt_1d[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(w, wpnt_1d.getWeight(), std::numeric_limits<double>::epsilon());
}

TEST(MathLib, WeightedPoint2D)
{
    std::array<double, 2> pnt;
    pnt[0] = 0.1;
    pnt[1] = 0.2;
    double w = 200.0;
    MathLib::WeightedPoint2D wpnt_2d(pnt, w);

    ASSERT_NEAR(pnt[0], wpnt_2d[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(pnt[1], wpnt_2d[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(w, wpnt_2d.getWeight(), std::numeric_limits<double>::epsilon());
}

TEST(MathLib, WeightedPoint3D)
{
    std::array<double, 3> pnt;
    pnt[0] = 0.1;
    pnt[1] = 0.2;
    pnt[2] = 0.3;
    double w = 300.0;
    MathLib::WeightedPoint3D wpnt_3d(pnt, w);

    ASSERT_NEAR(pnt[0], wpnt_3d[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(pnt[1], wpnt_3d[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(pnt[2], wpnt_3d[2], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(w, wpnt_3d.getWeight(), std::numeric_limits<double>::epsilon());
}
