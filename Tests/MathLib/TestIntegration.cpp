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

#include <vector>
#include <algorithm>

#include "MathLib/Integration/GaussLegendre.h"

namespace
{
double f1(double x) {return x*x;}
} //namespace

TEST(MathLib, IntegrationGaussLegendre)
{
    const double exact = 2./3.;

    ASSERT_EQ(0.0, MathLib::GaussLegendre::integrate(f1, 1));
    ASSERT_NEAR(exact, MathLib::GaussLegendre::integrate(f1, 2), exact*1e-5);
    ASSERT_NEAR(exact, MathLib::GaussLegendre::integrate(f1, 3), exact*1e-5);
    ASSERT_NEAR(exact, MathLib::GaussLegendre::integrate(f1, 4), exact*1e-5);
}


