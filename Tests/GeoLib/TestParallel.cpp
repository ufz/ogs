/**
 * \file
 * \date Feb 26, 2014
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "GeoLib/AnalyticalGeometry.h"

TEST(GeoLib, TestParallel)
{
    // parallel vectors
    Eigen::Vector3d v(0.0, 1.0, 2.0);
    Eigen::Vector3d w(0.0, 2.0, 4.0);
    EXPECT_TRUE(GeoLib::parallel(v, w));
    EXPECT_TRUE(GeoLib::parallel(w, v));

    v[1] = 0.0;
    w[1] = 0.0;
    EXPECT_TRUE(GeoLib::parallel(v, w));
    EXPECT_TRUE(GeoLib::parallel(w, v));

    // degenerate cases
    v[2] = 0.0;
    w[2] = 0.0;
    EXPECT_FALSE(GeoLib::parallel(v, w));
    EXPECT_FALSE(GeoLib::parallel(w, v));

    w[2] = 0.1;
    EXPECT_FALSE(GeoLib::parallel(v, w));
    EXPECT_FALSE(GeoLib::parallel(w, v));

    // non-parallel case
    v[1] = 0.1;
    EXPECT_FALSE(GeoLib::parallel(v, w));
    EXPECT_FALSE(GeoLib::parallel(w, v));

    // parallel vectors, opposite sense of direction
    v[0] = 0.0;
    v[1] = 1.0;
    v[2] = 2.0;
    w[0] = 0.0;
    w[1] = -2.0;
    w[2] = -4.0;
    EXPECT_TRUE(GeoLib::parallel(v, w));
    EXPECT_TRUE(GeoLib::parallel(w, v));
}
