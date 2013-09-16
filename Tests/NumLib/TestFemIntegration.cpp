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

#include <limits>

#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"

using namespace NumLib;

TEST(NumLib, FemIntegrationGaussRegular)
{
    const std::size_t integrationOrder = 2;
    const double eps = std::numeric_limits<double>::epsilon();

    // check position indices
    // dim = 1
    ASSERT_EQ((std::array<std::size_t, 1>({0u})), IntegrationGaussRegular<1>::getPositionIndices(integrationOrder, 0));
    ASSERT_EQ((std::array<std::size_t, 1>({1u})), IntegrationGaussRegular<1>::getPositionIndices(integrationOrder, 1));
    // dim = 2
    ASSERT_EQ((std::array<std::size_t, 2>({0u, 0u})), IntegrationGaussRegular<2>::getPositionIndices(integrationOrder, 0));
    ASSERT_EQ((std::array<std::size_t, 2>({0u, 1u})), IntegrationGaussRegular<2>::getPositionIndices(integrationOrder, 1));
    ASSERT_EQ((std::array<std::size_t, 2>({1u, 0u})), IntegrationGaussRegular<2>::getPositionIndices(integrationOrder, 2));
    ASSERT_EQ((std::array<std::size_t, 2>({1u, 1u})), IntegrationGaussRegular<2>::getPositionIndices(integrationOrder, 3));
    // dim = 3
    ASSERT_EQ((std::array<std::size_t, 3>({0u, 0u, 0u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 0));
    ASSERT_EQ((std::array<std::size_t, 3>({0u, 0u, 1u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 1));
    ASSERT_EQ((std::array<std::size_t, 3>({0u, 1u, 0u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 2));
    ASSERT_EQ((std::array<std::size_t, 3>({0u, 1u, 1u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 3));
    ASSERT_EQ((std::array<std::size_t, 3>({1u, 0u, 0u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 4));
    ASSERT_EQ((std::array<std::size_t, 3>({1u, 0u, 1u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 5));
    ASSERT_EQ((std::array<std::size_t, 3>({1u, 1u, 0u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 6));
    ASSERT_EQ((std::array<std::size_t, 3>({1u, 1u, 1u})), IntegrationGaussRegular<3>::getPositionIndices(integrationOrder, 7));

	// check coordinates
	// dim = 1
	// weight
	ASSERT_NEAR(1.0,
			IntegrationGaussRegular<1>::getWeightedPoint(integrationOrder, 0).getWeight(), eps);
	// pos
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<1>::getWeightedPoint(integrationOrder, 0)[0], eps);
	// weight
	ASSERT_NEAR(1.0,
			IntegrationGaussRegular<1>::getWeightedPoint(integrationOrder, 1).getWeight(), eps);
	// pos
	ASSERT_NEAR(-0.577350269189626,
			IntegrationGaussRegular<1>::getWeightedPoint(integrationOrder, 1)[0], eps);

	// dim = 2
	ASSERT_NEAR(1.0,
			IntegrationGaussRegular<2>::getWeightedPoint(integrationOrder, 0).getWeight(), eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<2>::getWeightedPoint(integrationOrder, 0)[0], eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<2>::getWeightedPoint(integrationOrder, 0)[1], eps);
	ASSERT_NEAR(1.0,
			IntegrationGaussRegular<2>::getWeightedPoint(integrationOrder, 1).getWeight(), eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<2>::getWeightedPoint(integrationOrder, 1)[0], eps);
	ASSERT_NEAR(-0.577350269189626,
			IntegrationGaussRegular<2>::getWeightedPoint(integrationOrder, 1)[1], eps);

	// dim = 3
	ASSERT_NEAR(1.0,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 0).getWeight(), eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 0)[0], eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 0)[1], eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 0)[2], eps);
	ASSERT_NEAR(1.0,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 1).getWeight(), eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 1)[0], eps);
	ASSERT_NEAR(0.577350269189626,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 1)[1], eps);
	ASSERT_NEAR(-0.577350269189626,
			IntegrationGaussRegular<3>::getWeightedPoint(integrationOrder, 1)[2], eps);

	// check other member functions
	IntegrationGaussRegular<1> q1;
	ASSERT_EQ(2u, q1.getIntegrationOrder());
	ASSERT_EQ(2u, q1.getNPoints());
	q1.setIntegrationOrder(3u);
	ASSERT_EQ(3u, q1.getIntegrationOrder());
	ASSERT_EQ(3u, q1.getNPoints());
	IntegrationGaussRegular<2> q2;
	ASSERT_EQ(2u, q2.getIntegrationOrder());
	ASSERT_EQ(4u, q2.getNPoints());
	IntegrationGaussRegular<3> q3;
	ASSERT_EQ(2u, q3.getIntegrationOrder());
	ASSERT_EQ(8u, q3.getNPoints());
}


