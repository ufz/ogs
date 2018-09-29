/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>
#include <memory>

#include <gtest/gtest.h>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "NumLib/Function/LinearInterpolationAlongPolyline.h"
#include "NumLib/Function/LinearInterpolationOnSurface.h"
#include "NumLib/Function/SpatialFunctionLinear.h"

#include "Tests/TestTools.h"

TEST(NumLib, SpatialFunctionLinear)
{
    // f(x,y,z) = 1 + 2x + 3y + 4z
    std::array<double,4> f_coeff = {{1, 2, 3, 4}};
    MathLib::LinearFunction<double,3> linear_f(f_coeff);

    NumLib::SpatialFunctionLinear f(linear_f);

    ASSERT_DOUBLE_EQ(1., f(GeoLib::Point(0,0,0)));
    ASSERT_DOUBLE_EQ(10., f(GeoLib::Point(1,1,1)));
    ASSERT_DOUBLE_EQ(-8, f(GeoLib::Point(-1,-1,-1)));
    for (std::size_t k(0); k < 5; ++k)
    {
        GeoLib::Point pt(0, 0, 0);
        for (unsigned i = 0; i < 3; ++i)
            pt[i] = (double) rand() - (RAND_MAX / 2.0);
        double expected = 1+2*pt[0]+3*pt[1]+4*pt[2];
        ASSERT_DOUBLE_EQ(expected, f(pt));
    }
}

TEST(NumLib, SpatialFunctionInterpolationPolyline)
{
    // create geometry
    GeoLib::Point pt1(0.0, 0.0, 0.0);
    GeoLib::Point pt2(10.0, 0.0, 0.0);
    std::vector<GeoLib::Point*> pnts = {&pt1, &pt2};
    GeoLib::Polyline ply0(pnts);
    ply0.addPoint(0);
    ply0.addPoint(1);

    // define a function
    const std::vector<std::size_t> vec_point_ids = {{0, 1}};
    const std::vector<double> vec_point_values = {{0., 100.}};
    NumLib::LinearInterpolationAlongPolyline interpolate(
            ply0, vec_point_ids, vec_point_values, std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::max());

    // normal
    for (unsigned k=0; k<10; k++)
        ASSERT_DOUBLE_EQ(k*10., interpolate(GeoLib::Point(k,0,0)));

    // failure
    // x
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(-1,0,0)));
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(11,0,0)));
    // y
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(0,1,0)));
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(0,-1,0)));
    // z
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(0,0,1)));
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(0,0,-1)));
}

TEST(NumLib, SpatialFunctionInterpolationSurface)
{
    // create geometry
    GeoLib::Point pt1(0.0, 0.0, 0.0);
    GeoLib::Point pt2(10.0, 0.0, 0.0);
    GeoLib::Point pt3(10.0, 10.0, 0.0);
    GeoLib::Point pt4(0.0, 10.0, 0.0);
    std::vector<GeoLib::Point*> pnts = {&pt1, &pt2, &pt3, &pt4};
    GeoLib::Surface sfc1(pnts);
    sfc1.addTriangle(0, 1, 2);
    sfc1.addTriangle(0, 2, 3);

    // define a function
    const std::vector<std::size_t> vec_point_ids = {{0, 1, 2, 3}};
    const std::vector<double> vec_point_values = {{0., 100., 100., 0.}};
    NumLib::LinearInterpolationOnSurface interpolate(
        sfc1, vec_point_ids, vec_point_values,
        std::numeric_limits<double>::max());

    // normal
    for (unsigned k=0; k<10; k++) {
        ASSERT_DOUBLE_EQ(k*10., interpolate(GeoLib::Point(k,0,0)));
        ASSERT_DOUBLE_EQ(k*10., interpolate(GeoLib::Point(k,5,0)));
        ASSERT_DOUBLE_EQ(k*10., interpolate(GeoLib::Point(k,10,0)));
    }

    // failure
    // x, y
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(-1,-1,0)));
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(11,-1,0)));
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(11,11,0)));
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(-1,11,0)));
    // z
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(0,0,1)));
    ASSERT_DOUBLE_EQ(std::numeric_limits<double>::max(), interpolate(GeoLib::Point(0,0,-1)));
}

