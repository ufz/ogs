/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <map>
#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "ProcessLib/TH2M/PhaseTransitionModels/PhaseTransitionNone.h"
#include "Tests/MaterialLib/TestMPL.h"
#include "Tests/TestTools.h"

TEST(ProcessLib, TH2MPhaseTransitionNone)
{
    std::stringstream m;

    auto const density_air = 1.2;
    auto const density_water = 999.5;

    auto const molar_mass_air = 0.02897;

    auto const specific_heat_capacity_air = 733.;
    auto const specific_heat_capacity_water = 4182.;

    auto const diffusion_coefficient_vapour = 0.0;

    auto const viscosity_air = 1.e-5;
    auto const viscosity_water = 1.e-3;

    m << "<medium>\n";
    m << "  <phases>\n";

    // gas phase
    m << "<phase>\n";
    m << "<type>Gas</type>\n";

    // gas phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("viscosity", viscosity_air);
    m << Tests::makeConstantPropertyElement("density", density_air);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_air);
    m << Tests::makeConstantPropertyElement("molar_mass", molar_mass_air);

    m << "</properties>\n";
    m << "</phase>\n";

    // liquid phase
    m << "<phase>\n";
    m << "<type>AqueousLiquid</type>\n";

    // liquid phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("viscosity", viscosity_water);
    m << Tests::makeConstantPropertyElement("density", density_water);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_water);
    m << "</properties>\n";
    m << "</phase>\n";
    m << "</phases> </medium>\n";

    std::shared_ptr<MaterialPropertyLib::Medium> const& medium =
        Tests::createTestMaterial(m.str());

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media{
        {0, medium}};

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    auto ptm = std::make_unique<ProcessLib::TH2M::PhaseTransitionNone>(media);

    double const pGR = 1000000.;
    double const pCap = 1000000.;
    double const T = 333.;

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)] = pGR;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = T;

    ptm->computeConstitutiveVariables(medium.get(), variable_array, pos, time,
                                      dt);
    auto const& cv = ptm->cv;

    // reference values
    double const rhoCGR = density_air;
    double const rhoWGR = 0.0;
    double const rhoCLR = 0.0;
    double const rhoWLR = density_water;
    double const xmCG = 1.0;
    double const xmWG = 0.0;
    double const dxmWG_dpGR = 0.0;
    double const dxmCG_dpGR = 0.0;
    double const dxmWG_dT = 0.0;
    double const dxmCG_dT = 0.0;
    double const hCG = 0.;
    double const hWG = 0.;
    double const hG = specific_heat_capacity_air * T;
    double const hL = specific_heat_capacity_water * T;
    double const uG = hG - pGR / density_air;
    double const uL = hL;
    double const muGR = viscosity_air;
    double const muLR = viscosity_water;

    ASSERT_NEAR(density_air, cv.rhoGR, 1.0e-10);
    ASSERT_NEAR(density_water, cv.rhoLR, 1.0e-10);
    ASSERT_NEAR(rhoCGR, cv.rhoCGR, 1.0e-10);
    ASSERT_NEAR(rhoWGR, cv.rhoWGR, 1.0e-10);
    ASSERT_NEAR(rhoCLR, cv.rhoCLR, 1.0e-10);
    ASSERT_NEAR(rhoWLR, cv.rhoWLR, 1.0e-10);
    ASSERT_NEAR(xmCG, cv.xmCG, 1.e-10);
    ASSERT_NEAR(xmWG, cv.xmWG, 1.e-10);
    ASSERT_NEAR(dxmWG_dpGR, cv.dxmWG_dpGR, 1.0e-10);
    ASSERT_NEAR(dxmCG_dpGR, cv.dxmCG_dpGR, 1.0e-10);
    ASSERT_NEAR(dxmWG_dT, cv.dxmWG_dT, 1.0e-10);
    ASSERT_NEAR(dxmCG_dT, cv.dxmCG_dT, 1.0e-10);
    ASSERT_NEAR(hCG, cv.hCG, 1.0e-9);
    ASSERT_NEAR(hWG, cv.hWG, 1.0e-9);
    ASSERT_NEAR(hG, cv.hG, 1.0e-10);
    ASSERT_NEAR(hL, cv.hL, 1.0e-10);
    ASSERT_NEAR(uG, cv.uG, 1.0e-10);
    ASSERT_NEAR(uL, cv.uL, 1.0e-10);
    ASSERT_NEAR(diffusion_coefficient_vapour, cv.diffusion_coefficient_vapour,
                1.0e-10);
    ASSERT_NEAR(muGR, cv.muGR, 1.0e-10);
    ASSERT_NEAR(muLR, cv.muLR, 1.0e-10);
}
