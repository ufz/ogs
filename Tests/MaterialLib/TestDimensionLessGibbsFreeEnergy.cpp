/**
 *  \copyright
 *   Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestDimensionLessGibbsFreeEnergy.cpp
 *
 * Created on December 9, 2016, 12:31 PM
 */
#include <gtest/gtest.h>

#include <memory>

#include "Tests/TestTools.h"

#include "MaterialLib/Fluid/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"

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
