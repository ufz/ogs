/**
 * \author Norihiro Watanabe
 * \date   2013-08-29
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

    EXPECT_EQ(0.0, MathLib::WeightedSum<MathLib::GaussLegendre<1>>::add(square));
    EXPECT_NEAR(2./3, MathLib::WeightedSum<MathLib::GaussLegendre<2>>::add(square),
            eps);
    EXPECT_NEAR(2./3, MathLib::WeightedSum<MathLib::GaussLegendre<3>>::add(square),
            eps);
    EXPECT_NEAR(2./3, MathLib::WeightedSum<MathLib::GaussLegendre<4>>::add(square),
            eps);

    auto const& cube = [](double const x) { return x * x * x; };
    EXPECT_NEAR(0.0, MathLib::WeightedSum<MathLib::GaussLegendre<1>>::add(cube),
                eps);
    EXPECT_NEAR(0.0, MathLib::WeightedSum<MathLib::GaussLegendre<2>>::add(cube),
                eps);
    EXPECT_NEAR(0.0, MathLib::WeightedSum<MathLib::GaussLegendre<3>>::add(cube),
                eps);
    EXPECT_NEAR(0.0, MathLib::WeightedSum<MathLib::GaussLegendre<4>>::add(cube),
                eps);
}


