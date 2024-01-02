/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/FormattingUtils.h"

TEST(MathLib, FormattingUtilsEigenDenseTypes)
{
    // integer 2-component column vector
    {
        Eigen::Vector2i v{5, 8};
        auto const actual = fmt::format("{}", v);
        std::string const expected = "5\n8";
        EXPECT_EQ(expected, actual);
    }

    // const double 3-component row vector
    {
        Eigen::RowVector3d const v{5.25, 8.5, -9.125};
        auto const actual = fmt::format("{}", v);
        std::string const expected = "  5.25    8.5 -9.125";
        EXPECT_EQ(expected, actual);
    }

    // 2 x 2 column major matrix
    {
        Eigen::Matrix2d m;
        m << 1, 2, 3, 4;
        auto const actual = fmt::format("{}", m);
        std::string const expected = "1 2\n3 4";
        EXPECT_EQ(expected, actual);
    }

    // 3 x 2 dynamic size row major matrix
    {
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m(
            3, 2);
        m << 11, 12, 21, 22, 31, 32;

        auto const actual = fmt::format("{}", m);
        std::string const expected = "11 12\n21 22\n31 32";
        EXPECT_EQ(expected, actual);
    }

    // array
    {
        Eigen::Array2d a{7, 13};

        auto const actual = fmt::format("{}", a);
        std::string const expected = " 7\n13";
        EXPECT_EQ(expected, actual);
    }

    // map
    {
        std::vector<double> d{3, 8};
        Eigen::Map<Eigen::Vector2d> v{d.data()};

        auto const actual = fmt::format("{}", v);
        std::string const expected = "3\n8";
        EXPECT_EQ(expected, actual);
    }
}

TEST(MathLib, FormattingUtilsEigenDenseExpressions)
{
    // mathematical expressions
    {
        Eigen::Array2i v{5, 8};
        auto const actual = fmt::format("{}", v * v);
        std::string const expected = "25\n64";
        EXPECT_EQ(expected, actual);
    }

    // reshape
    {
        auto const actual = fmt::format(
            "{}",
            Eigen::Vector4d::LinSpaced(1, 4).reshaped<Eigen::RowMajor>(2, 2));
        std::string const expected = "1 2\n3 4";
        EXPECT_EQ(expected, actual);
    }

    // block
    {
        Eigen::Matrix4d m = Eigen::Vector<double, 16>::LinSpaced(1, 16)
                                .reshaped<Eigen::ColMajor>(4, 4);

        auto const actual = fmt::format("{}", m.block<2, 3>(1, 1));
        std::string const expected = " 6 10 14\n 7 11 15";
        EXPECT_EQ(expected, actual);
    }
}
