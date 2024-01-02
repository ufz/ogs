/**
 * \file
 *  \copyright
 *   Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on December 9, 2016, 12:31 PM
 */
#include <gtest/gtest.h>

#include <array>
#include <memory>

#include "MaterialLib/Fluid/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"
#include "MaterialLib/Fluid/GibbsFreeEnergy/DimensionlessGibbsFreeEnergyRegion2.h"
#include "Tests/TestTools.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

TEST(MaterialGibbsFreeEnergy, checkDimensionLessGibbsFreeEnergyRegion1)
{
    DimensionLessGibbsFreeEnergyRegion1 gfeR1;
    const double T = 473.15;
    const double p = 4.e+7;
    const double T_ref = 1386;
    const double p_ref = 1.653e+7;

    const double tau = T_ref / T;
    const double pi = p / p_ref;

    ASSERT_NEAR(-0.946316961309, gfeR1.get_gamma(tau, pi), 1.e-10);
    ASSERT_NEAR(0.0849626110088, gfeR1.get_dgamma_dpi(tau, pi), 1.e-10);
    ASSERT_NEAR(1.360259836187, gfeR1.get_dgamma_dtau(tau, pi), 1.e-10);
    ASSERT_NEAR(-9.726737889528e-4, gfeR1.get_dgamma_dpi_dpi(tau, pi), 1.e-10);
    ASSERT_NEAR(-1.089765627732, gfeR1.get_dgamma_dtau_dtau(tau, pi), 1.e-10);
    ASSERT_NEAR(0.013153993854, gfeR1.get_dgamma_dtau_dpi(tau, pi), 1.e-10);
}

TEST(MaterialGibbsFreeEnergy, checkDimensionlessGibbsFreeEnergyRegion2)
{
    std::array const pressures = {1., 1., 5.e6, 20.e6, 99.e6};
    std::array const temperatures = {274.15, 623.15, 453.15, 624.15, 1072.15};

    std::array const expected_Gamma = {-6.487342246274503, -18.639988598929357,
                                       -813.70817666392134, -2.3869259565036165,
                                       -5.5887782055821402};
    std::array const expected_dGamma_dTau = {
        10.044285690295167, 12.75051476112089, -68588.395895402442,
        -5.543347937022638, 14.907400321522008};
    std::array const expected_dGamma_dPi = {
        999999.18256963324, 999999.98255043128, -2995.0842839942193,
        -0.094429398511927651, 0.0088374644316908044};

    for (std::size_t i = 0; i < pressures.size(); i++)
    {
        const double T_ref = 540;
        const double p_ref = 1.e+6;

        const double tau = T_ref / temperatures[i];
        const double pi = pressures[i] / p_ref;

        ASSERT_NEAR(expected_Gamma[i],
                    DimensionlessGibbsFreeEnergyRegion2::getGamma(tau, pi),
                    1.e-10);
        ASSERT_NEAR(expected_dGamma_dTau[i],
                    DimensionlessGibbsFreeEnergyRegion2::getdGammadTau(tau, pi),
                    1.e-10);
        ASSERT_NEAR(expected_dGamma_dPi[i],
                    DimensionlessGibbsFreeEnergyRegion2::getdGammadPi(tau, pi),
                    1.e-10);
    }
}
