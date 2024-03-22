/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <map>
#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "ProcessLib/TH2M/ConstitutiveRelations/NoPhaseTransition.h"
#include "Tests/MaterialLib/TestMPL.h"
#include "Tests/TestTools.h"

TEST(ProcessLib, TH2MNoPhaseTransition)
{
    using namespace ProcessLib::TH2M;
    using namespace ProcessLib::TH2M::ConstitutiveRelations;

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

    // solid phase
    m << "<phase>\n";
    m << "<type>Solid</type>\n";
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("density", 2e3);
    m << "</properties>\n";
    m << "</phase>\n";

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
    m << "</phases>\n";
    m << "<properties>\n";
    m << "  <property>\n";
    m << "      <name>tortuosity</name>\n";
    m << "      <type>Constant</type>\n";
    m << "      <value>1</value>\n";
    m << "  </property>\n";
    m << "</properties> </medium>\n";

    std::shared_ptr<MaterialPropertyLib::Medium> const& medium =
        Tests::createTestMaterial(m.str());
    ProcessLib::TH2M::MediaData media_data{*medium};

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media{
        {0, medium}};

    MaterialPropertyLib::VariableArray variable_array;
    ProcessLib::ConstitutiveRelations::SpaceTimeData x_t{
        {},
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};

    auto ptm = std::make_unique<NoPhaseTransition>(media);

    double const pGR = 1000000.;
    double const pCap = 1000000.;
    double const T = 333.;

    variable_array.gas_phase_pressure = pGR;
    variable_array.capillary_pressure = pCap;
    variable_array.temperature = T;

    ProcessLib::TH2M::ConstitutiveRelations::PureLiquidDensityData rhoWLR;
    ProcessLib::TH2M::ConstitutiveRelations::PureLiquidDensityModel
        rhoWLR_model;
    rhoWLR_model.eval(x_t, media_data, GasPressureData{pGR},
                      CapillaryPressureData{pGR}, TemperatureData{T, T},
                      rhoWLR);
    ASSERT_NEAR(density_water, rhoWLR(), 1e-10);

    ProcessLib::TH2M::ConstitutiveRelations::ViscosityData viscosity;
    ProcessLib::TH2M::ConstitutiveRelations::EnthalpyData enthalpy;
    PhaseTransitionData cv;
    ptm->eval(x_t, media_data, GasPressureData{pGR}, CapillaryPressureData{pGR},
              TemperatureData{T, T}, rhoWLR, viscosity, enthalpy, cv);

    // reference values
    double const rhoCGR = density_air;
    double const rhoWGR = 0.0;
    double const rhoCLR = 0.0;
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
    ASSERT_NEAR(xmCG, 1. - cv.xmWG, 1.e-10);
    ASSERT_NEAR(xmWG, cv.xmWG, 1.e-10);
    ASSERT_NEAR(dxmWG_dpGR, cv.dxmWG_dpGR, 1.0e-10);
    ASSERT_NEAR(dxmCG_dpGR, -cv.dxmWG_dpGR, 1.0e-10);
    ASSERT_NEAR(dxmWG_dT, cv.dxmWG_dT, 1.0e-10);
    ASSERT_NEAR(dxmCG_dT, -cv.dxmWG_dT, 1.0e-10);
    ASSERT_NEAR(hCG, cv.hCG, 1.0e-9);
    ASSERT_NEAR(hWG, cv.hWG, 1.0e-9);
    ASSERT_NEAR(hG, enthalpy.h_G, 1.0e-10);
    ASSERT_NEAR(hL, enthalpy.h_L, 1.0e-10);
    ASSERT_NEAR(uG, cv.uG, 1.0e-10);
    ASSERT_NEAR(uL, cv.uL, 1.0e-10);
    ASSERT_NEAR(diffusion_coefficient_vapour, cv.diffusion_coefficient_vapour,
                1.0e-10);
    ASSERT_NEAR(muGR, viscosity.mu_GR, 1.0e-10);
    ASSERT_NEAR(muLR, viscosity.mu_LR, 1.0e-10);
}
