/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <limits>
#include <type_traits>
#include <gtest/gtest.h>
#include <logog/include/logog.hpp>
#include "MathLib/Nonlinear/Root1D.h"

double f(double x)
{
    return x*x-1;
}

template <typename T>
class MathLibRegulaFalsi : public ::testing::Test
{
};

namespace NL = MathLib::Nonlinear;
using RegulaFalsiTypes = ::testing::Types<NL::Unmodified, NL::Illinois, NL::Pegasus, NL::AndersonBjorck>;

TYPED_TEST_CASE(MathLibRegulaFalsi, RegulaFalsiTypes);

TYPED_TEST(MathLibRegulaFalsi, QuadraticFunction)
{
    auto rf = NL::makeRegulaFalsi<TypeParam>(f, -0.1, 1.1);
    double old_range = rf.getRange();

    DBUG(" 0 -- x ~ %23.16g, range = %23.16g", rf.getResult(), old_range);

    for (unsigned n=0; n<10; ++n)
    {
        rf.step(1);
        double range = rf.getRange();
        // expect that the interval of the root search shrinks
        EXPECT_GT(old_range, range);
        old_range = range;
        DBUG("%2i -- x ~ %23.16g, range = %23.16g", n+1, rf.getResult(), range);

        if (range < std::numeric_limits<double>::epsilon()) break;
    }

    auto const error = std::abs(f(rf.getResult()));

    if (!std::is_same<NL::Unmodified, TypeParam>::value) {
        EXPECT_GT(std::numeric_limits<double>::epsilon(), old_range);
        EXPECT_GT(std::numeric_limits<double>::epsilon(), error);
    } else {
        // The unmodified regula falsi method converges very slowly.
        EXPECT_GT(100.0*std::numeric_limits<double>::epsilon(), error);
    }
}
