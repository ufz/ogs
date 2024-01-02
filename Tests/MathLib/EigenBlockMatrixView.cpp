/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MathLib/EigenBlockMatrixView.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

TEST(MathLib_EigenBlockMatriView, 1by1)
{
    Eigen::Matrix<double, 4, 4> expected =
        Eigen::Vector4d::Constant(42).asDiagonal();

    // fixed size
    Eigen::Matrix<double, 1, 1> A{42};
    ASSERT_EQ(expected, MathLib::eigenBlockMatrixView<4>(A));

    // dynamic size
    Eigen::MatrixXd B(1, 1);
    B << 42;
    ASSERT_TRUE(expected == MathLib::eigenBlockMatrixView<4>(B));
}

TEST(MathLib_EigenBlockMatriView, 1byN)
{
    Eigen::Matrix<double, 6, 3> expected;
    // clang-format off
    expected << 1, 0, 0,
                2, 0, 0,
                0, 1, 0,
                0, 2, 0,
                0, 0, 1,
                0, 0, 2;
    // clang-format on

    Eigen::Vector2d v{1, 2};
    ASSERT_TRUE(expected == MathLib::eigenBlockMatrixView<3>(v));

    // dynamic
    Eigen::MatrixXd w(2, 1);
    w << 1, 2;
    ASSERT_TRUE(expected == MathLib::eigenBlockMatrixView<3>(w));
}

TEST(MathLib_EigenBlockMatriView, Nby1)
{
    Eigen::Matrix<double, 3, 6> expected;
    // clang-format off
    expected << 1, 2, 0, 0, 0, 0,
                0, 0, 1, 2, 0, 0,
                0, 0, 0, 0, 1, 2;
    // clang-format on

    Eigen::RowVector2d v{1, 2};
    ASSERT_TRUE(expected == MathLib::eigenBlockMatrixView<3>(v));

    // dynamic
    Eigen::MatrixXd w(1, 2);
    w << 1, 2;
    ASSERT_TRUE(expected == MathLib::eigenBlockMatrixView<3>(w));
}

TEST(MathLib_EigenBlockMatriView, NbyM)
{
    Eigen::Matrix<double, 6, 6> expected;
    // clang-format off
    expected << 1, 2, 0, 0, 0, 0,
                3, 4, 0, 0, 0, 0,
                0, 0, 1, 2, 0, 0,
                0, 0, 3, 4, 0, 0,
                0, 0, 0, 0, 1, 2,
                0, 0, 0, 0, 3, 4;
    // clang-format on

    Eigen::Matrix2d A;
    A << 1, 2, 3, 4;
    ASSERT_TRUE(expected == MathLib::eigenBlockMatrixView<3>(A));

    // dynamic
    Eigen::MatrixXd B(2, 2);
    B << 1, 2, 3, 4;
    ASSERT_TRUE(expected == MathLib::eigenBlockMatrixView<3>(B));
}
