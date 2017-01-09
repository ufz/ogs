/**
 * @file TestComputeRotationMatrix.cpp
 * @date 2015-04-23
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include "GeoLib/AnalyticalGeometry.h"

auto test3equal = [](double a, double b, double c, double const* result)
{
    EXPECT_EQ(a, result[0]);
    EXPECT_EQ(b, result[1]);
    EXPECT_EQ(c, result[2]);

    delete[] result;
};

TEST(GeoLib, ComputeRotationMatrixToXYnegative)
{
    MathLib::Vector3 const n(0.0, -1.0, 0.0);
    MathLib::DenseMatrix<double> rot_mat(3,3,0.0);

    GeoLib::computeRotationMatrixToXY(n, rot_mat);
    EXPECT_EQ(1.0, rot_mat(0,0));
    EXPECT_EQ(0.0, rot_mat(0,1));
    EXPECT_EQ(0.0, rot_mat(0,2));
    EXPECT_EQ(0.0, rot_mat(1,0));
    EXPECT_EQ(0.0, rot_mat(1,1));
    EXPECT_EQ(1.0, rot_mat(1,2));
    EXPECT_EQ(0.0, rot_mat(2,0));
    EXPECT_EQ(-1.0, rot_mat(2,1));
    EXPECT_EQ(0.0, rot_mat(2,2));

    MathLib::Vector3 const x(0.0,1.0,0.0);
    test3equal(0, 0, -1, rot_mat*x.getCoords());

    MathLib::Vector3 const x0(10.0,1.0,0.0);
    test3equal(10, 0, -1, rot_mat*x0.getCoords());

    MathLib::Vector3 const x1(10.0,0.0,10.0);
    test3equal(10, 10, 0, rot_mat*x1.getCoords());
}

TEST(GeoLib, ComputeRotationMatrixToXYpositive)
{
    MathLib::Vector3 const n(0.0, 1.0, 0.0);
    MathLib::DenseMatrix<double> rot_mat(3,3,0.0);

    GeoLib::computeRotationMatrixToXY(n, rot_mat);
    EXPECT_EQ(1.0, rot_mat(0,0));
    EXPECT_EQ(0.0, rot_mat(0,1));
    EXPECT_EQ(0.0, rot_mat(0,2));
    EXPECT_EQ(0.0, rot_mat(1,0));
    EXPECT_EQ(0.0, rot_mat(1,1));
    EXPECT_EQ(-1.0, rot_mat(1,2));
    EXPECT_EQ(0.0, rot_mat(2,0));
    EXPECT_EQ(1.0, rot_mat(2,1));
    EXPECT_EQ(0.0, rot_mat(2,2));

    MathLib::Vector3 const x(0.0,1.0,0.0);
    test3equal(0, 0, 1, rot_mat*x.getCoords());

    MathLib::Vector3 const x0(10.0,1.0,0.0);
    test3equal(10, 0, 1, rot_mat*x0.getCoords());

    MathLib::Vector3 const x1(10.0,0.0,10.0);
    test3equal(10, -10, 0, rot_mat*x1.getCoords());
}

