// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "MathLib/Integration/GaussLegendre.h"

namespace
{
double square(double const x)
{
    return x * x;
}
}  // namespace

TEST(MathLib, IntegrationGaussLegendre)
{
    double const eps = 10 * std::numeric_limits<double>::epsilon();

    EXPECT_EQ(0.0,
              MathLib::WeightedSum<MathLib::GaussLegendre<1>>::add(square));
    EXPECT_NEAR(2. / 3,
                MathLib::WeightedSum<MathLib::GaussLegendre<2>>::add(square),
                eps);
    EXPECT_NEAR(2. / 3,
                MathLib::WeightedSum<MathLib::GaussLegendre<3>>::add(square),
                eps);
    EXPECT_NEAR(2. / 3,
                MathLib::WeightedSum<MathLib::GaussLegendre<4>>::add(square),
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
