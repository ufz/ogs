// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <map>
#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "ProcessLib/TH2M/ConstitutiveRelations/NoPhaseTransition.h"
#include "Tests/MaterialLib/TestMPL.h"
#include "Tests/TestTools.h"

TEST(ProcessLib, TH2MNoPhaseTransition)
{
    namespace CR = ProcessLib::TH2M::ConstitutiveRelations;

    std::stringstream m;

    auto const density_air = 1.2;
    auto const density_water = 999.5;

    auto const molar_mass_air = 0.02897;

    auto const specific_heat_capacity_air = 733.;
    auto const specific_heat_capacity_water = 4182.;

    auto const diffusion_coefficient_vapour = 0.0;

    m << "<medium>\n";
    m << "  <phases>\n";

    // solid phase
    m << "<phase>\n";
    m << "<type>Solid</type>\n";
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("density", 2e3);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity", 0);
    m << Tests::makeConstantPropertyElement("thermal_expansivity", 0);
    m << "</properties>\n";
    m << "</phase>\n";

    // gas phase
    m << "<phase>\n";
    m << "<type>Gas</type>\n";

    // gas phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("density", density_air);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_air);
    m << Tests::makeConstantPropertyElement("molar_mass", molar_mass_air);
    m << Tests::makeConstantPropertyElement("viscosity", 0);

    m << "</properties>\n";
    m << "</phase>\n";

    // liquid phase
    m << "<phase>\n";
    m << "<type>AqueousLiquid</type>\n";

    // liquid phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("density", density_water);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_water);
    m << Tests::makeConstantPropertyElement("viscosity", 0);
    m << "</properties>\n";
    m << "</phase>\n";
    m << "</phases>\n";
    m << "<properties>\n";
    m << "  <property>\n";
    m << "      <name>tortuosity</name>\n";
    m << "      <type>Constant</type>\n";
    m << "      <value>1</value>\n";
    m << "  </property>\n";
    m << Tests::makeConstantPropertyElement("saturation", 0);
    m << Tests::makeConstantPropertyElement("permeability", 0);
    m << Tests::makeConstantPropertyElement("relative_permeability", 0);
    m << Tests::makeConstantPropertyElement(
        "relative_permeability_nonwetting_phase", 0);
    m << Tests::makeConstantPropertyElement("bishops_effective_stress", 0);
    m << Tests::makeConstantPropertyElement("porosity", 0);
    m << Tests::makeConstantPropertyElement("biot_coefficient", 0);
    m << Tests::makeConstantPropertyElement("thermal_conductivity", 0);
    m << "</properties> </medium>\n";

    std::shared_ptr<MaterialPropertyLib::Medium> const& medium =
        Tests::createTestMaterial(m.str());
    CR::MediaData media_data{*medium};

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media{
        {0, medium}};

    MaterialPropertyLib::VariableArray variable_array;
    CR::SpaceTimeData x_t{{},
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};

    auto ptm = std::make_unique<CR::NoPhaseTransition>(media);

    double const pGR = 1000000.;
    double const pCap = 1000000.;
    double const T = 333.;

    variable_array.gas_phase_pressure = pGR;
    variable_array.capillary_pressure = pCap;
    variable_array.temperature = T;

    CR::PureLiquidDensityData rhoWLR;
    CR::PureLiquidDensityModel rhoWLR_model;
    rhoWLR_model.eval(x_t, media_data, CR::GasPressureData{pGR},
                      CR::CapillaryPressureData{pGR}, CR::TemperatureData{T, T},
                      rhoWLR);
    ASSERT_NEAR(density_water, rhoWLR(), 1e-10);

    CR::FluidEnthalpyData enthalpy;
    CR::MassMoleFractionsData mass_mole_fractions;
    CR::FluidDensityData fluid_density;
    CR::VapourPartialPressureData vapour_pressure;
    CR::ConstituentDensityData constituent_density;
    CR::PhaseTransitionData cv;
    ptm->eval(x_t, media_data, CR::GasPressureData{pGR},
              CR::CapillaryPressureData{pGR}, CR::TemperatureData{T, T}, rhoWLR,
              enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
              constituent_density, cv);

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

    ASSERT_NEAR(density_air, fluid_density.rho_GR, 1.0e-10);
    ASSERT_NEAR(density_water, fluid_density.rho_LR, 1.0e-10);
    ASSERT_NEAR(rhoCGR, constituent_density.rho_C_GR, 1.0e-10);
    ASSERT_NEAR(rhoWGR, constituent_density.rho_W_GR, 1.0e-10);
    ASSERT_NEAR(rhoCLR, constituent_density.rho_C_LR, 1.0e-10);
    ASSERT_NEAR(xmCG, mass_mole_fractions.xmCG, 1.e-10);
    ASSERT_NEAR(xmWG, 1 - mass_mole_fractions.xmCG, 1.e-10);
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
}
