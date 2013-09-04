/**
 * @file WeightedPoint.cpp
 * @author Thomas Fischer
 * @date Sep 4, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// google test
#include "gtest/gtest.h"

#include "TemplateWeightedPoint.h"

TEST(MathLib, WeightedPoint1D)
{
	double pnt = 0.5;
	double w = 100.0;
	MathLib::WeightedPoint1D wpnt_1d({pnt}, w);

	ASSERT_NEAR(wpnt_1d[0], pnt, std::numeric_limits<double>::min());
	ASSERT_NEAR(wpnt_1d.getWeight(), w, std::numeric_limits<double>::min());
}

TEST(MathLib, WeightedPoint2D)
{
	double pnt[2] = {0.1, 0.2};
	double w = 200.0;
	MathLib::WeightedPoint2D wpnt_2d({pnt[0], pnt[1]}, w);

	ASSERT_NEAR(wpnt_2d[0], pnt[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(wpnt_2d[1], pnt[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(wpnt_2d.getWeight(), w, std::numeric_limits<double>::min());
}

TEST(MathLib, WeightedPoint3D)
{
	double pnt[3] = {0.1, 0.2, 0.3};
	double w = 300.0;
	MathLib::WeightedPoint3D wpnt_3d({pnt[0], pnt[1], pnt[2]}, w);

	ASSERT_NEAR(wpnt_3d[0], pnt[0], std::numeric_limits<double>::min());
	ASSERT_NEAR(wpnt_3d[1], pnt[1], std::numeric_limits<double>::min());
	ASSERT_NEAR(wpnt_3d[2], pnt[2], std::numeric_limits<double>::min());
	ASSERT_NEAR(wpnt_3d.getWeight(), w, std::numeric_limits<double>::min());
}
