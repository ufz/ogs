/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "MaterialLib/MPL/Property.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, PropertyFromVector1Elem)
{
    std::vector const values{3.14};
    auto const prop = MPL::fromVector(values);

    using V = double;
    ASSERT_TRUE(std::holds_alternative<V>(prop));
    ASSERT_EQ(3.14, get<V>(prop));
}

TEST(MaterialPropertyLib, PropertyFromVector2Elems)
{
    std::vector const values{2.71, 3.14};
    auto const prop = MPL::fromVector(values);

    using V = Eigen::Vector2d;
    ASSERT_TRUE(std::holds_alternative<V>(prop));
    ASSERT_EQ((V{2.71, 3.14}), get<V>(prop));
}

TEST(MaterialPropertyLib, PropertyFromVector3Elems)
{
    std::vector const values{2.71, 3.14, 1.602e-19};
    auto const prop = MPL::fromVector(values);

    using V = Eigen::Vector3d;
    ASSERT_TRUE(std::holds_alternative<V>(prop));
    ASSERT_EQ((V{2.71, 3.14, 1.602e-19}), get<V>(prop));
}

TEST(MaterialPropertyLib, PropertyFromVector4Elems)
{
    std::vector const values{2.71, 3.14,  //
                             1.602e-19, 6.626e-34};
    auto const prop = MPL::fromVector(values);

    using V = Eigen::Matrix2d;
    V expected;
    expected(0, 0) = 2.71;
    expected(0, 1) = 3.14;  // row major storage order of the values vector
    expected(1, 0) = 1.602e-19;
    expected(1, 1) = 6.626e-34;

    ASSERT_TRUE(std::holds_alternative<V>(prop));
    ASSERT_EQ(expected, get<V>(prop));
}

TEST(MaterialPropertyLib, PropertyFromVector6Elems)
{
    std::vector const values{2.71,      3.14,   1.602e-19,
                             6.626e-34, 2.99e8, 1. / 137.};
    auto const prop = MPL::fromVector(values);

    using V = Eigen::Vector<double, 6>;
    V expected;
    expected << 2.71, 3.14, 1.602e-19, 6.626e-34, 2.99e8, 1. / 137.;

    ASSERT_TRUE(std::holds_alternative<V>(prop));
    ASSERT_EQ(expected, get<V>(prop));
}

TEST(MaterialPropertyLib, PropertyFromVector9Elems)
{
    std::vector const values{2.71,      3.14,     1.602e-19,  //
                             6.626e-34, 2.99e8,   1. / 137.,  //
                             1.38e-23,  9.11e-31, 2.002};
    auto const prop = MPL::fromVector(values);

    using V = Eigen::Matrix3d;
    V expected;
    expected.row(0) << 2.71, 3.14, 1.602e-19;
    expected.row(1) << 6.626e-34, 2.99e8, 1. / 137.;
    expected.row(2) << 1.38e-23, 9.11e-31, 2.002;

    ASSERT_TRUE(std::holds_alternative<V>(prop));
    ASSERT_EQ(expected, get<V>(prop));
}

TEST(MaterialPropertyLib, PropertyFromVectorInvalidElemCount)
{
    std::vector const values{2.71,      3.14,     1.602e-19,  //
                             6.626e-34, 2.99e8,   1. / 137.,  //
                             1.38e-23,  9.11e-31, 2.002};

    for (std::size_t invalid_size : {0, 5, 7, 8})
    {
        std::vector const vs(values.begin(), values.begin() + invalid_size);
        EXPECT_ANY_THROW(MPL::fromVector(vs))
            << "exception expected for invalid input size of " << invalid_size;
    }
}
