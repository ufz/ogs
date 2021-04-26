/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <map>
#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "ProcessLib/TH2M/PhaseTransitionModels/PhaseTransitionEvaporation.h"
#include "Tests/MaterialLib/TestMPL.h"
#include "Tests/TestTools.h"

TEST(ProcessLib, TH2MPhaseTransitionEvaporation)
{
    std::stringstream m;

    auto const density_air = 1.2;
    auto const density_water = 999.5;

    auto const molar_mass_air = 0.02897;
    auto const molar_mass_water = 0.018053;

    auto const specific_heat_capacity_air = 733.;
    auto const specific_heat_capacity_water = 4182.;

    auto const vapour_pressure_water = 18900.;
    auto const specific_latent_heat_water = 2258000.;

    auto const diffusion_coefficient_vapour = 2.6e-6;

    auto const thermal_conductivity_air = 0.2;
    auto const thermal_conductivity_water = 0.6;
    auto const viscosity_air = 1.e-5;
    auto const viscosity_water = 1.e-3;

    m << "<medium>\n";
    m << "  <phases>\n";

    // gas phase
    m << "<phase>\n";
    m << "<type>Gas</type>\n";
    m << " <components>\n";

    // gas phase component Air
    m << "  <component>\n";
    m << "   <name>Air</name>\n";
    m << "   <properties>\n";
    m << Tests::makeConstantPropertyElement("molar_mass", molar_mass_air);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_air);
    m << "   </properties>\n";
    m << "  </component>\n";

    // gas phase component Water Vapour
    m << "  <component>\n";
    m << "   <name>Vapour</name>\n";
    m << "   <properties>\n";
    m << Tests::makeConstantPropertyElement("molar_mass", molar_mass_water);
    m << Tests::makeConstantPropertyElement("vapour_pressure",
                                            vapour_pressure_water);
    m << Tests::makeConstantPropertyElement("diffusion",
                                            diffusion_coefficient_vapour);
    m << Tests::makeConstantPropertyElement("specific_latent_heat",
                                            specific_latent_heat_water);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_water);
    m << "   </properties>\n";
    m << "  </component>\n";
    m << " </components>\n";

    // gas phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("viscosity", viscosity_air);
    m << Tests::makeConstantPropertyElement("density", density_air);
    m << Tests::makeConstantPropertyElement("thermal_conductivity",
                                            thermal_conductivity_air);
    m << "</properties>\n";
    m << "</phase>\n";

    // liquid phase
    m << "<phase>\n";
    m << "<type>AqueousLiquid</type>\n";

    // liquid phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("viscosity", viscosity_water);
    m << Tests::makeConstantPropertyElement("density", density_water);
    m << Tests::makeConstantPropertyElement("thermal_conductivity",
                                            thermal_conductivity_water);
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

    std::unique_ptr<ProcessLib::TH2M::PhaseTransitionModels> ptm =
        std::make_unique<ProcessLib::TH2M::PhaseTransitionEvaporation>(media);

    double const pGR = 1000000.;
    double const pCap = 1000000.;
    double const T = 333.;

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)] = pGR;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = T;

    ptm->getConstitutiveVariables(medium.get(), variable_array, pos, time, dt);

    // reference values
    double const rhoCGR = 1.185858524394711;
    double const rhoWGR = 0.01414147560528907;
    double const xmCG = 0.9882154369955924;
    double const xmWG = 0.01178456300440756;
    double const dxmWG_dpGR = 0.;
    double const dxmCG_dpGR = -0.;
    double const dxmWG_dT = -0.0003043909214558457;
    double const dxmCG_dT = 0.0003043909214558457;

    double const hCG = specific_heat_capacity_air * T;
    double const hWG =
        specific_heat_capacity_water * T + specific_latent_heat_water;

    double const hG = xmCG * hCG + xmWG * (hWG);
    double const hL = specific_heat_capacity_water * T;
    double const uG = hG - pGR / density_air;
    double const uL = hL;
    double const muGR = viscosity_air;
    double const lambdaGR = thermal_conductivity_air;
    double const muLR = viscosity_water;
    double const lambdaLR = thermal_conductivity_water;

    ASSERT_NEAR(density_air, ptm->rhoGR, 1.0e-10);
    ASSERT_NEAR(density_water, ptm->rhoLR, 1.0e-10);
    ASSERT_NEAR(rhoCGR, ptm->rhoCGR, 1.0e-10);
    ASSERT_NEAR(rhoWGR, ptm->rhoWGR, 1.0e-10);
    ASSERT_NEAR(xmCG, ptm->xmCG, 1.e-10);
    ASSERT_NEAR(xmWG, ptm->xmWG, 1.e-10);
    ASSERT_NEAR(dxmWG_dpGR, ptm->dxmWG_dpGR, 1.0e-10);
    ASSERT_NEAR(dxmCG_dpGR, ptm->dxmCG_dpGR, 1.0e-10);
    ASSERT_NEAR(dxmWG_dT, ptm->dxmWG_dT, 1.0e-10);
    ASSERT_NEAR(dxmCG_dT, ptm->dxmCG_dT, 1.0e-10);
    ASSERT_NEAR(hCG, ptm->hCG, 1.0e-09);
    ASSERT_NEAR(hWG, ptm->hWG, 1.0e-09);
    ASSERT_NEAR(hG, ptm->hG, 1.0e-09);
    ASSERT_NEAR(hL, ptm->hL, 1.0e-09);
    ASSERT_NEAR(uG, ptm->uG, 1.0e-09);
    ASSERT_NEAR(uL, ptm->uL, 1.0e-09);
    ASSERT_NEAR(diffusion_coefficient_vapour, ptm->diffusion_coefficient_vapour,
                1.0e-10);
    ASSERT_NEAR(muGR, ptm->muGR, 1.0e-10);
    ASSERT_NEAR(lambdaGR, ptm->lambdaGR, 1.0e-10);
    ASSERT_NEAR(muLR, ptm->muLR, 1.0e-10);
    ASSERT_NEAR(lambdaLR, ptm->lambdaLR, 1.0e-10);
}