/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <numeric>
#include <random>

#include "NumLib/DOF/LocalDOF.h"

static double testLocalDOFRandom()
{
    std::random_device device;
    std::default_random_engine engine(device());
    std::uniform_real_distribution<double> distribution(100, 200);
    return distribution(engine);
}

template <unsigned NPoints>
struct MockShapeFunction
{
    static constexpr unsigned NPOINTS = NPoints;
};

TEST(NumLib_LocalDOF, Scalar)
{
    using N1 = MockShapeFunction<2>;
    using N2 = MockShapeFunction<3>;
    using N3 = MockShapeFunction<4>;

    auto const value = testLocalDOFRandom();
    auto const x = Eigen::VectorXd::LinSpaced(9, value, value + 8).eval();

    auto const [x1, x2, x3] = NumLib::localDOF<N1, N2, N3>(x);

#ifndef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
    static_assert(decltype(x1)::RowsAtCompileTime == 2);
    static_assert(decltype(x2)::RowsAtCompileTime == 3);
    static_assert(decltype(x3)::RowsAtCompileTime == 4);
#endif

    ASSERT_EQ(x.rows(), x1.rows() + x2.rows() + x3.rows());
    ASSERT_EQ(2, x1.rows());
    ASSERT_EQ(3, x2.rows());
    ASSERT_EQ(4, x3.rows());

    EXPECT_DOUBLE_EQ(x[0], x1[0]);
    EXPECT_DOUBLE_EQ(x[1], x1[1]);
    EXPECT_DOUBLE_EQ(x[2], x2[0]);
    EXPECT_DOUBLE_EQ(x[3], x2[1]);
    EXPECT_DOUBLE_EQ(x[4], x2[2]);
    EXPECT_DOUBLE_EQ(x[5], x3[0]);
    EXPECT_DOUBLE_EQ(x[6], x3[1]);
    EXPECT_DOUBLE_EQ(x[7], x3[2]);
    EXPECT_DOUBLE_EQ(x[8], x3[3]);
}

TEST(NumLib_LocalDOF, Vectorial)
{
    using N1 = MockShapeFunction<2>;
    using N2 = MockShapeFunction<3>;
    using N3 = MockShapeFunction<4>;

    constexpr auto num_dof_total = 2 + 3 * 2 + 4;

    // also works with std::vector
    auto const value = testLocalDOFRandom();
    std::vector<double> x(num_dof_total);
    std::iota(x.begin(), x.end(), value);

    {
        auto const [x1, x2, x3] =
            NumLib::localDOF<N1, NumLib::Vectorial<N2, 2>, N3>(x);

#ifndef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
        static_assert(decltype(x1)::RowsAtCompileTime == 2);
        static_assert(decltype(x2)::RowsAtCompileTime == 3 * 2);
        static_assert(decltype(x3)::RowsAtCompileTime == 4);
#endif

        ASSERT_EQ(num_dof_total, x1.rows() + x2.rows() + x3.rows());
        ASSERT_EQ(2, x1.rows());
        ASSERT_EQ(3 * 2, x2.rows());
        ASSERT_EQ(4, x3.rows());

        EXPECT_DOUBLE_EQ(x[0], x1[0]);
        EXPECT_DOUBLE_EQ(x[1], x1[1]);

        EXPECT_DOUBLE_EQ(x[2], x2[0]);
        EXPECT_DOUBLE_EQ(x[3], x2[1]);
        EXPECT_DOUBLE_EQ(x[4], x2[2]);
        EXPECT_DOUBLE_EQ(x[5], x2[3]);
        EXPECT_DOUBLE_EQ(x[6], x2[4]);
        EXPECT_DOUBLE_EQ(x[7], x2[5]);

        EXPECT_DOUBLE_EQ(x[8], x3[0]);
        EXPECT_DOUBLE_EQ(x[9], x3[1]);
        EXPECT_DOUBLE_EQ(x[10], x3[2]);
        EXPECT_DOUBLE_EQ(x[11], x3[3]);
    }

    // same but with vectorial d.o.f. first
    {
        auto const [x1, x2, x3] =
            NumLib::localDOF<NumLib::Vectorial<N2, 2>, N1, N3>(x);

#ifndef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
        static_assert(decltype(x1)::RowsAtCompileTime == 3 * 2);
        static_assert(decltype(x2)::RowsAtCompileTime == 2);
        static_assert(decltype(x3)::RowsAtCompileTime == 4);
#endif

        ASSERT_EQ(num_dof_total, x1.rows() + x2.rows() + x3.rows());
        ASSERT_EQ(3 * 2, x1.rows());
        ASSERT_EQ(2, x2.rows());
        ASSERT_EQ(4, x3.rows());

        EXPECT_DOUBLE_EQ(x[0], x1[0]);
        EXPECT_DOUBLE_EQ(x[1], x1[1]);
        EXPECT_DOUBLE_EQ(x[2], x1[2]);
        EXPECT_DOUBLE_EQ(x[3], x1[3]);
        EXPECT_DOUBLE_EQ(x[4], x1[4]);
        EXPECT_DOUBLE_EQ(x[5], x1[5]);

        EXPECT_DOUBLE_EQ(x[6], x2[0]);
        EXPECT_DOUBLE_EQ(x[7], x2[1]);

        EXPECT_DOUBLE_EQ(x[8], x3[0]);
        EXPECT_DOUBLE_EQ(x[9], x3[1]);
        EXPECT_DOUBLE_EQ(x[10], x3[2]);
        EXPECT_DOUBLE_EQ(x[11], x3[3]);
    }
}
