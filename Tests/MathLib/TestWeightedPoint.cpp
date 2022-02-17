/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "MathLib/WeightedPoint.h"

TEST(MathLib, WeightedPoint0D)
{
    std::array<double, 0> const pnt{};
    double const w = 50.0;
    MathLib::WeightedPoint const wpnt_0d(pnt, w);

    ASSERT_EQ(0, wpnt_0d.getDimension());
    ASSERT_EQ(w, wpnt_0d.getWeight());
}

TEST(MathLib, WeightedPoint1D)
{
    std::array<double, 1> pnt{0.5};
    double const w = 100.0;
    MathLib::WeightedPoint const wpnt_1d(pnt, w);

    ASSERT_EQ(1, wpnt_1d.getDimension());
    ASSERT_EQ(pnt[0], wpnt_1d[0]);
    ASSERT_EQ(w, wpnt_1d.getWeight());
}

TEST(MathLib, WeightedPoint2D)
{
    std::array<double, 2> const pnt{0.1, 0.2};
    double const w = 200.0;
    MathLib::WeightedPoint const wpnt_2d(pnt, w);

    ASSERT_EQ(2, wpnt_2d.getDimension());
    ASSERT_EQ(pnt[0], wpnt_2d[0]);
    ASSERT_EQ(pnt[1], wpnt_2d[1]);
    ASSERT_EQ(w, wpnt_2d.getWeight());
}

TEST(MathLib, WeightedPoint3D)
{
    std::array<double, 3> const pnt{0.1, 0.2, 0.3};
    double const w = 300.0;
    MathLib::WeightedPoint const wpnt_3d(pnt, w);

    ASSERT_EQ(3, wpnt_3d.getDimension());
    ASSERT_EQ(pnt[0], wpnt_3d[0]);
    ASSERT_EQ(pnt[1], wpnt_3d[1]);
    ASSERT_EQ(pnt[2], wpnt_3d[2]);
    ASSERT_EQ(w, wpnt_3d.getWeight());
}

TEST(MathLib, WeightedPointEquality3D)
{
    MathLib::WeightedPoint const orig{std::array{0.1, 0.2, 0.3}, 0.4};

    {
        MathLib::WeightedPoint const same{std::array{0.1, 0.2, 0.3}, 0.4};
        ASSERT_EQ(orig, same);
    }

    {
        MathLib::WeightedPoint const diff_x{std::array{0.15, 0.2, 0.3}, 0.4};
        ASSERT_NE(orig, diff_x);
    }

    {
        MathLib::WeightedPoint const diff_y{std::array{0.1, 0.25, 0.3}, 0.4};
        ASSERT_NE(orig, diff_y);
    }

    {
        MathLib::WeightedPoint const diff_z{std::array{0.1, 0.2, 0.35}, 0.4};
        ASSERT_NE(orig, diff_z);
    }

    {
        MathLib::WeightedPoint const diff_w{std::array{0.1, 0.2, 0.3}, 0.45};
        ASSERT_NE(orig, diff_w);
    }

    {
        MathLib::WeightedPoint const diff_dim0{std::array<double, 0>{}, 0.4};
        ASSERT_NE(orig, diff_dim0);
    }

    {
        MathLib::WeightedPoint const diff_dim1{std::array{0.1}, 0.4};
        ASSERT_NE(orig, diff_dim1);
    }

    {
        MathLib::WeightedPoint const diff_dim2{std::array{0.1, 0.2}, 0.4};
        ASSERT_NE(orig, diff_dim2);
    }
}

TEST(MathLib, WeightedPointEquality2D)
{
    MathLib::WeightedPoint const orig{std::array{0.1, 0.2}, 0.4};

    {
        MathLib::WeightedPoint const same{std::array{0.1, 0.2}, 0.4};
        ASSERT_EQ(orig, same);
    }

    {
        MathLib::WeightedPoint const diff_x{std::array{0.15, 0.2}, 0.4};
        ASSERT_NE(orig, diff_x);
    }

    {
        MathLib::WeightedPoint const diff_y{std::array{0.1, 0.25}, 0.4};
        ASSERT_NE(orig, diff_y);
    }

    {
        MathLib::WeightedPoint const diff_w{std::array{0.1, 0.2}, 0.45};
        ASSERT_NE(orig, diff_w);
    }

    {
        MathLib::WeightedPoint const diff_dim0{std::array<double, 0>{}, 0.4};
        ASSERT_NE(orig, diff_dim0);
    }

    {
        MathLib::WeightedPoint const diff_dim1{std::array{0.1}, 0.4};
        ASSERT_NE(orig, diff_dim1);
    }
}

TEST(MathLib, WeightedPointEquality1D)
{
    MathLib::WeightedPoint const orig{std::array{0.1}, 0.4};

    {
        MathLib::WeightedPoint const same{std::array{0.1}, 0.4};
        ASSERT_EQ(orig, same);
    }

    {
        MathLib::WeightedPoint const diff_x{std::array{0.15}, 0.4};
        ASSERT_NE(orig, diff_x);
    }

    {
        MathLib::WeightedPoint const diff_w{std::array{0.1}, 0.45};
        ASSERT_NE(orig, diff_w);
    }

    {
        MathLib::WeightedPoint const diff_dim0{std::array<double, 0>{}, 0.4};
        ASSERT_NE(orig, diff_dim0);
    }
}

TEST(MathLib, WeightedPointEquality0D)
{
    MathLib::WeightedPoint const orig{0.4};

    {
        MathLib::WeightedPoint const same{0.4};
        ASSERT_EQ(orig, same);
    }

    {
        MathLib::WeightedPoint const diff_w{0.45};
        ASSERT_NE(orig, diff_w);
    }
}
