/*!
   \file  TestDensity.cpp
   \brief Test classes of density.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include "MaterialLib/Solid/Density/LinearSolidDensityModel.h"
#include "MaterialLib/Fluid/Density/IdealGasLaw.h"
#include "MaterialLib/ConstantDensityModel.h"
#include "MaterialLib/Density.h"

namespace
{

TEST(Material, checkDensity)
{
    //-- Solid --------------------------------------------------------------
    // Constant
    constexpr double rho = 2080.;
    MaterialLib::Density<MaterialLib::ConstantDensityModel> s_density(rho);
    ASSERT_NEAR(rho, s_density.getDensity(), 1.e-10);

    // Linear
    MaterialLib::Density<MaterialLib::LinearSolidDensityModel>
               lin_density(20., 1000., 100., 900.);
    const double lin_den_expected = (900. - 1000.)/(100. - 20.)
                                     * (50.-20.) + 1000.;
    ASSERT_NEAR(lin_den_expected, lin_density.getDensity(50.), 1.e-10);

    //-- Fluid --------------------------------------------------------------
    constexpr double molar_air = 28.96;
    MaterialLib::Density<MaterialLib::IdealGasLaw> air_density(molar_air);
    constexpr double T = 290.;
    constexpr double p = 1.e+5;
    constexpr double R = 8315.41;
    constexpr double expected_air_dens = molar_air * p /(R * T);
    ASSERT_NEAR(expected_air_dens, air_density.getDensity(T, p), 1.e-10);
}

}
