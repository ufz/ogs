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

#include <tuple>
#include <limits>

#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"

using namespace NumLib;

TEST(NumLib, FemIntegrationGaussRegular)
{
    const std::size_t nGaussLevel = 2;
    const double eps = std::numeric_limits<double>::epsilon();

    // check position indices
    // dim = 1
    ASSERT_EQ(std::make_tuple(0u, 0u, 0u), IntegrationGaussRegular<1>::getPosition(nGaussLevel, 0));
    ASSERT_EQ(std::make_tuple(1u, 0u, 0u), IntegrationGaussRegular<1>::getPosition(nGaussLevel, 1));
    // dim = 2
    ASSERT_EQ(std::make_tuple(0u, 0u, 0u), IntegrationGaussRegular<2>::getPosition(nGaussLevel, 0));
    ASSERT_EQ(std::make_tuple(0u, 1u, 0u), IntegrationGaussRegular<2>::getPosition(nGaussLevel, 1));
    ASSERT_EQ(std::make_tuple(1u, 0u, 0u), IntegrationGaussRegular<2>::getPosition(nGaussLevel, 2));
    ASSERT_EQ(std::make_tuple(1u, 1u, 0u), IntegrationGaussRegular<2>::getPosition(nGaussLevel, 3));
    // dim = 3
    ASSERT_EQ(std::make_tuple(0u, 0u, 0u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 0));
    ASSERT_EQ(std::make_tuple(0u, 0u, 1u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 1));
    ASSERT_EQ(std::make_tuple(0u, 1u, 0u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 2));
    ASSERT_EQ(std::make_tuple(0u, 1u, 1u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 3));
    ASSERT_EQ(std::make_tuple(1u, 0u, 0u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 4));
    ASSERT_EQ(std::make_tuple(1u, 0u, 1u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 5));
    ASSERT_EQ(std::make_tuple(1u, 1u, 0u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 6));
    ASSERT_EQ(std::make_tuple(1u, 1u, 1u), IntegrationGaussRegular<3>::getPosition(nGaussLevel, 7));

    // check coordinates
    double x[3];
    // dim = 1
    ASSERT_NEAR(1.0, IntegrationGaussRegular<1>::getPoint(nGaussLevel, 0, x), eps);
    ASSERT_NEAR(0.577350269189626, x[0], eps);
    ASSERT_NEAR(1.0, IntegrationGaussRegular<1>::getPoint(nGaussLevel, 1, x), eps);
    ASSERT_NEAR(-0.577350269189626, x[0], eps);
    // dim = 2
    ASSERT_NEAR(1.0, IntegrationGaussRegular<2>::getPoint(nGaussLevel, 0, x), eps);
    ASSERT_NEAR(0.577350269189626, x[0], eps);
    ASSERT_NEAR(0.577350269189626, x[1], eps);
    ASSERT_NEAR(1.0, IntegrationGaussRegular<2>::getPoint(nGaussLevel, 1, x), eps);
    ASSERT_NEAR(0.577350269189626, x[0], eps);
    ASSERT_NEAR(-0.577350269189626, x[1], eps);
    // dim = 3
    ASSERT_NEAR(1.0, IntegrationGaussRegular<3>::getPoint(nGaussLevel, 0, x), eps);
    ASSERT_NEAR(0.577350269189626, x[0], eps);
    ASSERT_NEAR(0.577350269189626, x[1], eps);
    ASSERT_NEAR(0.577350269189626, x[2], eps);
    ASSERT_NEAR(1.0, IntegrationGaussRegular<3>::getPoint(nGaussLevel, 1, x), eps);
    ASSERT_NEAR(0.577350269189626, x[0], eps);
    ASSERT_NEAR(0.577350269189626, x[1], eps);
    ASSERT_NEAR(-0.577350269189626, x[2], eps);

    // check other member functions
    IntegrationGaussRegular<1> q1;
    ASSERT_EQ(2u, q1.getSamplingLevel());
    ASSERT_EQ(2u, q1.getNPoints());
    q1.setSamplingLevel(3u);
    ASSERT_EQ(3u, q1.getSamplingLevel());
    ASSERT_EQ(3u, q1.getNPoints());
    IntegrationGaussRegular<2> q2;
    ASSERT_EQ(2u, q2.getSamplingLevel());
    ASSERT_EQ(4u, q2.getNPoints());
    IntegrationGaussRegular<3> q3;
    ASSERT_EQ(2u, q3.getSamplingLevel());
    ASSERT_EQ(8u, q3.getNPoints());
}


