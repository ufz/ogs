/**
 * @file TestCheckParallelism.cpp
 * @date Feb 26, 2014
 *
 * @copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "AnalyticalGeometry.h"

TEST(GeoLib, TestCheckParallelism)
{
	// parallel vectors
	double v[3] = {0.0, 1.0, 2.0};
	double w[3] = {0.0, 2.0, 4.0};
	EXPECT_TRUE(GeoLib::checkParallelism(v,w));

	v[1] = 0.0;
	w[1] = 0.0;
	EXPECT_TRUE(GeoLib::checkParallelism(v,w));

	// degenerate cases
	v[2] = 0.0;
	w[2] = 0.0;
	EXPECT_FALSE(GeoLib::checkParallelism(v,w));

	w[2] = 0.1;
	EXPECT_FALSE(GeoLib::checkParallelism(v,w));

	// non-parallel case
	v[1] = 0.1;
	EXPECT_FALSE(GeoLib::checkParallelism(v,w));

	// parallel vectors, opposite sense of direction
	v[0] = 0.0; v[1] = 1.0; v[2] = 2.0;
	w[0] = 0.0; w[1] = -2.0; w[2] = -4.0;
	EXPECT_TRUE(GeoLib::checkParallelism(v,w));
}

