/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-29
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <limits>

#include "NumLib/Fem/Integration/IntegrationGaussLegendreRegular.h"
#include "Tests/TestTools.h"

using namespace NumLib;

TEST(NumLib, FemIntegrationGaussLegendreRegular)
{
    const std::size_t integrationOrder = 2;
    const double eps = std::numeric_limits<double>::epsilon();

    // check position indices
    // dim = 1
    {
        std::size_t expected[1] = {0u};
        ASSERT_ARRAY_EQ(expected,
                        IntegrationGaussLegendreRegular<1>::getPositionIndices(
                            integrationOrder, 0),
                        1u);
        expected[0] = 1u;
        ASSERT_ARRAY_EQ(expected,
                        IntegrationGaussLegendreRegular<1>::getPositionIndices(
                            integrationOrder, 1),
                        1u);
    }
    // dim = 2
    {
        std::size_t expected[2] = {0u, 0u};
        ASSERT_ARRAY_EQ(expected,
                        IntegrationGaussLegendreRegular<2>::getPositionIndices(
                            integrationOrder, 0),
                        2);
        expected[1] = 1u;
        ASSERT_ARRAY_EQ(expected,
                        IntegrationGaussLegendreRegular<2>::getPositionIndices(
                            integrationOrder, 1),
                        2);
        expected[0] = 1u;
        expected[1] = 0u;
        ASSERT_ARRAY_EQ(expected,
                        IntegrationGaussLegendreRegular<2>::getPositionIndices(
                            integrationOrder, 2),
                        2);
        expected[0] = 1u;
        expected[1] = 1u;
        ASSERT_ARRAY_EQ(expected,
                        IntegrationGaussLegendreRegular<2>::getPositionIndices(
                            integrationOrder, 3),
                        2);
    }
    // dim = 3
    {
        std::size_t expected[3] = {0u, 0u, 0u};
        for (std::size_t i(0); i <= 1; i++)
        {
            expected[0] = i;
            for (std::size_t j(0); j <= 1; j++)
            {
                expected[1] = j;
                for (std::size_t k(0); k <= 1; k++)
                {
                    expected[2] = k;
                    const std::size_t l(i * 4 + j * 2 + k);
                    ASSERT_ARRAY_EQ(
                        expected,
                        IntegrationGaussLegendreRegular<3>::getPositionIndices(
                            integrationOrder, l),
                        3);
                }
            }
        }
    }
    // check coordinates
    // dim = 1
    // weight
    ASSERT_NEAR(1.0,
                IntegrationGaussLegendreRegular<1>::getWeightedPoint(
                    integrationOrder, 0)
                    .getWeight(),
                eps);
    // pos
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<1>::getWeightedPoint(
                    integrationOrder, 0)[0],
                eps);
    // weight
    ASSERT_NEAR(1.0,
                IntegrationGaussLegendreRegular<1>::getWeightedPoint(
                    integrationOrder, 1)
                    .getWeight(),
                eps);
    // pos
    ASSERT_NEAR(-0.577350269189626,
                IntegrationGaussLegendreRegular<1>::getWeightedPoint(
                    integrationOrder, 1)[0],
                eps);

    // dim = 2
    ASSERT_NEAR(1.0,
                IntegrationGaussLegendreRegular<2>::getWeightedPoint(
                    integrationOrder, 0)
                    .getWeight(),
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<2>::getWeightedPoint(
                    integrationOrder, 0)[0],
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<2>::getWeightedPoint(
                    integrationOrder, 0)[1],
                eps);
    ASSERT_NEAR(1.0,
                IntegrationGaussLegendreRegular<2>::getWeightedPoint(
                    integrationOrder, 1)
                    .getWeight(),
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<2>::getWeightedPoint(
                    integrationOrder, 1)[0],
                eps);
    ASSERT_NEAR(-0.577350269189626,
                IntegrationGaussLegendreRegular<2>::getWeightedPoint(
                    integrationOrder, 1)[1],
                eps);

    // dim = 3
    ASSERT_NEAR(1.0,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 0)
                    .getWeight(),
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 0)[0],
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 0)[1],
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 0)[2],
                eps);
    ASSERT_NEAR(1.0,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 1)
                    .getWeight(),
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 1)[0],
                eps);
    ASSERT_NEAR(0.577350269189626,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 1)[1],
                eps);
    ASSERT_NEAR(-0.577350269189626,
                IntegrationGaussLegendreRegular<3>::getWeightedPoint(
                    integrationOrder, 1)[2],
                eps);

    // check other member functions
    IntegrationGaussLegendreRegular<1> q1;
    ASSERT_EQ(2u, q1.getIntegrationOrder());
    ASSERT_EQ(2u, q1.getNumberOfPoints());
    q1.setIntegrationOrder(3u);
    ASSERT_EQ(3u, q1.getIntegrationOrder());
    ASSERT_EQ(3u, q1.getNumberOfPoints());
    IntegrationGaussLegendreRegular<2> q2;
    ASSERT_EQ(2u, q2.getIntegrationOrder());
    ASSERT_EQ(4u, q2.getNumberOfPoints());
    IntegrationGaussLegendreRegular<3> q3;
    ASSERT_EQ(2u, q3.getIntegrationOrder());
    ASSERT_EQ(8u, q3.getNumberOfPoints());
}
