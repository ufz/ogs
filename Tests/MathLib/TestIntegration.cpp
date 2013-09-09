/**
 * \author Norihiro Watanabe
 * \date   2013-08-29
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "MathLib/Integration/GaussLegendre.h"

namespace
{
    double square(double const x)
    {
        return x*x;
    }
}

TEST(MathLib, IntegrationGaussLegendre)
{
    double const eps = 10 * std::numeric_limits<double>::epsilon();

    EXPECT_EQ(0.0, MathLib::GaussLegendre::integrate(square, 1));
    EXPECT_NEAR(2./3, MathLib::GaussLegendre::integrate(square, 2),
            eps);
    EXPECT_NEAR(2./3, MathLib::GaussLegendre::integrate(square, 3),
            eps);
    EXPECT_NEAR(2./3, MathLib::GaussLegendre::integrate(square, 4),
            eps);

    auto const& cube = [](double const x){ return x*x*x; };
    EXPECT_EQ(0.0, MathLib::GaussLegendre::integrate(cube, 1));
    EXPECT_EQ(0.0, MathLib::GaussLegendre::integrate(cube, 2));
    EXPECT_EQ(0.0, MathLib::GaussLegendre::integrate(cube, 3));
    EXPECT_EQ(0.0, MathLib::GaussLegendre::integrate(cube, 4));
}


