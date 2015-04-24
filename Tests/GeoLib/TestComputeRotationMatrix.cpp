/**
 * @file TestComputeRotationMatrix.cpp
 * @date 2015-04-23
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"
#include "AnalyticalGeometry.h"

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
	MathLib::Vector3 const result(rot_mat*x.getCoords());
	EXPECT_EQ(0.0, result[0]);
	EXPECT_EQ(0.0, result[1]);
	EXPECT_EQ(-1.0, result[2]);

	MathLib::Vector3 const x0(10.0,1.0,0.0);
	MathLib::Vector3 const r0(rot_mat*x0.getCoords());
	EXPECT_EQ(10.0, r0[0]);
	EXPECT_EQ(0.0, r0[1]);
	EXPECT_EQ(-1.0, r0[2]);

	MathLib::Vector3 const x1(10.0,0.0,10.0);
	MathLib::Vector3 const r1(rot_mat*x1.getCoords());
	EXPECT_EQ(10.0, r1[0]);
	EXPECT_EQ(10.0, r1[1]);
	EXPECT_EQ(0.0, r1[2]);
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
	MathLib::Vector3 const result(rot_mat*x.getCoords());
	EXPECT_EQ(0.0, result[0]);
	EXPECT_EQ(0.0, result[1]);
	EXPECT_EQ(1.0, result[2]);

	MathLib::Vector3 const x0(10.0,1.0,0.0);
	MathLib::Vector3 const r0(rot_mat*x0.getCoords());
	EXPECT_EQ(10.0, r0[0]);
	EXPECT_EQ(0.0, r0[1]);
	EXPECT_EQ(1.0, r0[2]);

	MathLib::Vector3 const x1(10.0,0.0,10.0);
	MathLib::Vector3 const r1(rot_mat*x1.getCoords());
	EXPECT_EQ(10.0, r1[0]);
	EXPECT_EQ(-10.0, r1[1]);
	EXPECT_EQ(0.0, r1[2]);
}

