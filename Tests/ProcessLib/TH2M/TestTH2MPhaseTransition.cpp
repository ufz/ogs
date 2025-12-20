// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <map>
#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"
#include "ProcessLib/TH2M/ConstitutiveRelations/PhaseTransition.h"
#include "ProcessLib/TH2M/ConstitutiveRelations/PureLiquidDensity.h"
#include "Tests/MaterialLib/TestMPL.h"
#include "Tests/TestTools.h"

static const double molar_mass_air = 0.02897;
static const double molar_mass_water = 0.018053;
static const double specific_heat_capacity_air = 733.;
static const double specific_heat_capacity_water = 4182.;
static const double specific_latent_heat_water = 2258000.;
static const double diffusion_coefficient_vapour = 2.6e-6;
static const double thermal_conductivity_air = 0.2;
static const double thermal_conductivity_water = 0.6;

// In case of constant gas density
static const double constant_gas_density = 1.;

// Water density EOS
static const double rhoLR_ref = 1000.;
static const double T_ref = 288.15;
static const double pLR_ref = 1.01235e5;
static const double slope_T = -4.4e-4;
static const double slope_pLR = 4.65e-10;

std::string MediumDefinition(const bool density_is_constant)
{
    std::stringstream m;

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

    m << "   <property>\n";
    m << "     <name>vapour_pressure</name>\n";
    m << "     <type>ClausiusClapeyron</type>\n";
    m << "     <triple_temperature>273.16</triple_temperature>\n";
    m << "     <triple_pressure>836.21849</triple_pressure>\n";
    m << "     <critical_temperature>647.1</critical_temperature>\n";
    m << "     <critical_pressure>26016399.68391</critical_pressure>\n";
    m << "     <reference_temperature>373.15</reference_temperature>\n";
    m << "     <reference_pressure>101325</reference_pressure>\n";
    m << "   </property> \n";

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
    m << Tests::makeConstantPropertyElement("thermal_conductivity",
                                            thermal_conductivity_air);
    m << Tests::makeConstantPropertyElement("viscosity", 0);

    if (density_is_constant)
    {
        m << Tests::makeConstantPropertyElement("density",
                                                constant_gas_density);
    }
    else
    {
        m << "<property>\n";
        m << "<name>density</name>\n";
        m << "<type>IdealGasLawBinaryMixture</type>\n";
        m << "</property>\n";
    }

    m << "</properties>\n";
    m << "</phase>\n";

    // liquid phase
    m << "<phase>\n";
    m << "<type>AqueousLiquid</type>\n";

    // liquid phase components
    m << "<components>\n";
    m << "<component>\n";
    m << "<name>A</name>\n";
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("molar_mass", 0.);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity", 0.);
    m << Tests::makeConstantPropertyElement("henry_coefficient", 0.);
    m << Tests::makeConstantPropertyElement("diffusion", 0.);
    m << Tests::makeConstantPropertyElement("specific_latent_heat", 0.);
    m << "</properties>\n";
    m << "</component>\n";

    m << "<component>\n";
    m << "<name>W</name>\n";
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("molar_mass", 0.);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_water);
    m << "</properties>\n";
    m << "</component>\n";
    m << "</components>\n";

    // liquid phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_water);
    m << Tests::makeConstantPropertyElement("viscosity", 0);

    m << "  <property>\n";
    m << "      <name>density</name>\n";
    m << "      <type>Linear</type>\n";
    m << "      <reference_value>" << rhoLR_ref << "</reference_value>\n";
    m << "      <independent_variable>\n";
    m << "          <variable_name>temperature</variable_name>\n";
    m << "          <reference_condition>" << T_ref
      << "                    </reference_condition>\n";
    m << "          <slope>" << slope_T << "</slope>\n";
    m << "      </independent_variable>\n";
    m << "      <independent_variable>\n";
    m << "          <variable_name>liquid_phase_pressure</variable_name>\n";
    m << "          <reference_condition>" << pLR_ref
      << "                      </reference_condition>\n";
    m << "          <slope>" << slope_pLR << "</slope>\n";
    m << "      </independent_variable>\n";
    m << "  </property>\n";

    m << Tests::makeConstantPropertyElement("thermal_conductivity",
                                            thermal_conductivity_water);
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

    return m.str();
}

TEST(ProcessLib, TH2MPhaseTransition)
{
    namespace CR = ProcessLib::TH2M::ConstitutiveRelations;

    bool const density_is_constant = false;
    std::shared_ptr<MaterialPropertyLib::Medium> const& medium =
        Tests::createTestMaterial(MediumDefinition(density_is_constant));
    CR::MediaData media_data{*medium};

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media{
        {0, medium}};

    auto ptm = std::make_unique<CR::PhaseTransition>(media);

    auto const count = 200000;
    auto const pGR_min = 100000.;
    auto const pCap_min = 5000.;
    auto const T_min = 280.;

    auto const pGR_max = 200000.;
    auto const pCap_max = 5000000.;
    auto const T_max = 480.;

    auto const dpGR = (pGR_max - pGR_min) / count;
    auto const dpCap = (pCap_max - pCap_min) / count;
    auto const dT = (T_max - T_min) / count;
    auto const eps_pGR = pGR_max * 5.e-4;
    auto const eps_pCap = pCap_max * 2.e-5;
    auto const eps_T = T_max * 2.e-5;

    CR::PureLiquidDensityData rhoWLR;
    CR::PureLiquidDensityModel rhoWLR_model;

    CR::FluidEnthalpyData enthalpy;
    CR::MassMoleFractionsData mass_mole_fractions;
    CR::FluidDensityData fluid_density;
    CR::VapourPartialPressureData vapour_pressure;
    CR::ConstituentDensityData constituent_density;
    CR::PhaseTransitionData cv;
    CR::SpaceTimeData x_t{{},
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};

    for (std::size_t i = 0; i <= count; i++)
    {
        auto pGR = pGR_min + i * dpGR;
        auto pCap = pCap_min + i * dpCap;
        auto T = T_min + i * dT;
        CR::GasPressureData p_GR_data{pGR};
        CR::CapillaryPressureData p_cap_data{pCap};
        CR::TemperatureData T_data{T, T};

        // Perturb gas pressure
        rhoWLR_model.eval(x_t, media_data, CR::GasPressureData{pGR + eps_pGR},
                          p_cap_data, T_data, rhoWLR);
        ptm->eval(x_t, media_data, CR::GasPressureData{pGR + eps_pGR},
                  p_cap_data, T_data, rhoWLR, enthalpy, mass_mole_fractions,
                  fluid_density, vapour_pressure, constituent_density, cv);

        auto xmWG_plus = 1 - mass_mole_fractions.xmCG;
        auto rhoGR_plus = fluid_density.rho_GR;
        auto rhoCGR_plus = constituent_density.rho_C_GR;
        auto rhoWGR_plus = constituent_density.rho_W_GR;

        rhoWLR_model.eval(x_t, media_data, CR::GasPressureData{pGR - eps_pGR},
                          p_cap_data, T_data, rhoWLR);
        ptm->eval(x_t, media_data, CR::GasPressureData{pGR - eps_pGR},
                  p_cap_data, T_data, rhoWLR, enthalpy, mass_mole_fractions,
                  fluid_density, vapour_pressure, constituent_density, cv);

        auto xmWG_minus = 1 - mass_mole_fractions.xmCG;
        auto rhoGR_minus = fluid_density.rho_GR;
        auto rhoCGR_minus = constituent_density.rho_C_GR;
        auto rhoWGR_minus = constituent_density.rho_W_GR;

        // Unperturbed primary variables
        rhoWLR_model.eval(x_t, media_data, CR::GasPressureData{pGR}, p_cap_data,
                          T_data, rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        // Central difference derivatives
        auto const dxmWG_dpGR = (xmWG_plus - xmWG_minus) / (2. * eps_pGR);
        auto const drhoGR_dpGR = (rhoGR_plus - rhoGR_minus) / (2. * eps_pGR);
        auto const drhoCGR_dpGR = (rhoCGR_plus - rhoCGR_minus) / (2. * eps_pGR);
        auto const drhoWGR_dpGR = (rhoWGR_plus - rhoWGR_minus) / (2. * eps_pGR);

        // Test composition derivatives w.r.t. gas pressure
        auto const tolerance_dpGR = 1e-4;

        if ((1. - xmWG_plus) > tolerance_dpGR)
        {
            ASSERT_NEAR(dxmWG_dpGR, cv.dxmWG_dpGR, tolerance_dpGR);
            ASSERT_NEAR(drhoGR_dpGR, cv.drho_GR_dp_GR, tolerance_dpGR);
            ASSERT_NEAR(drhoCGR_dpGR, cv.drho_C_GR_dp_GR, tolerance_dpGR);
            ASSERT_NEAR(drhoWGR_dpGR, cv.drho_W_GR_dp_GR, tolerance_dpGR);
        }

        // Perturb capillary pressure
        rhoWLR_model.eval(x_t, media_data, p_GR_data,
                          CR::CapillaryPressureData{pCap + eps_pCap}, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data,
                  CR::CapillaryPressureData{pCap + eps_pCap}, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_plus = 1 - mass_mole_fractions.xmCG;
        rhoGR_plus = fluid_density.rho_GR;
        rhoCGR_plus = constituent_density.rho_C_GR;
        rhoWGR_plus = constituent_density.rho_W_GR;

        rhoWLR_model.eval(x_t, media_data, p_GR_data,
                          CR::CapillaryPressureData{pCap - eps_pCap}, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data,
                  CR::CapillaryPressureData{pCap - eps_pCap}, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_minus = 1 - mass_mole_fractions.xmCG;
        rhoGR_minus = fluid_density.rho_GR;
        rhoCGR_minus = constituent_density.rho_C_GR;
        rhoWGR_minus = constituent_density.rho_W_GR;

        // Unperturbed primary variables
        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        // Central difference derivatives
        auto const dxmWG_dpCap = (xmWG_plus - xmWG_minus) / (2. * eps_pCap);
        auto const drhoGR_dpCap = (rhoGR_plus - rhoGR_minus) / (2. * eps_pCap);
        auto const drhoCGR_dpCap =
            (rhoCGR_plus - rhoCGR_minus) / (2. * eps_pCap);
        auto const drhoWGR_dpCap =
            (rhoWGR_plus - rhoWGR_minus) / (2. * eps_pCap);

        // Test composition derivatives w.r.t. capillary pressure
        auto const tolerance_dpCap = 1e-10;

        if ((1. - xmWG_plus) > tolerance_dpCap)
        {
            ASSERT_NEAR(dxmWG_dpCap, cv.dxmWG_dpCap, tolerance_dpCap);
            ASSERT_NEAR(drhoGR_dpCap, cv.drho_GR_dp_cap, tolerance_dpCap);
            ASSERT_NEAR(drhoCGR_dpCap, cv.drho_C_GR_dp_cap, tolerance_dpCap);
            ASSERT_NEAR(drhoWGR_dpCap, cv.drho_W_GR_dp_cap, tolerance_dpCap);
        }

        // Perturb temperature
        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data,
                          CR::TemperatureData{T + eps_T, T}, rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data,
                  CR::TemperatureData{T + eps_T, T}, rhoWLR, enthalpy,
                  mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_plus = 1 - mass_mole_fractions.xmCG;
        rhoGR_plus = fluid_density.rho_GR;
        rhoCGR_plus = constituent_density.rho_C_GR;
        rhoWGR_plus = constituent_density.rho_W_GR;

        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data,
                          CR::TemperatureData{T - eps_T, T}, rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data,
                  CR::TemperatureData{T - eps_T, T}, rhoWLR, enthalpy,
                  mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_minus = 1 - mass_mole_fractions.xmCG;
        rhoGR_minus = fluid_density.rho_GR;
        rhoCGR_minus = constituent_density.rho_C_GR;
        rhoWGR_minus = constituent_density.rho_W_GR;

        // Unperturbed primary variables
        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        // Central difference derivatives
        auto const dxmWG_dT = (xmWG_plus - xmWG_minus) / (2. * eps_T);
        auto const drhoGR_dT = (rhoGR_plus - rhoGR_minus) / (2. * eps_T);
        auto const drhoCGR_dT = (rhoCGR_plus - rhoCGR_minus) / (2. * eps_T);
        auto const drhoWGR_dT = (rhoWGR_plus - rhoWGR_minus) / (2. * eps_T);

        // Test composition derivatives w.r.t. temperature
        auto const tolerance_dT = 1.e-4;

        if ((1. - xmWG_plus) > tolerance_dT)
        {
            ASSERT_NEAR(dxmWG_dT, cv.dxmWG_dT, tolerance_dT);
            ASSERT_NEAR(drhoGR_dT, cv.drho_GR_dT, tolerance_dT);
            ASSERT_NEAR(drhoCGR_dT, cv.drho_C_GR_dT, tolerance_dT);
            ASSERT_NEAR(drhoWGR_dT, cv.drho_W_GR_dT, tolerance_dT);
        }
        // Test mixture composition
        // Vapour mole fraction is the quotient of vapour pressure and gas
        // phase pressure
        auto const refence_xnWG =
            std::clamp(vapour_pressure.pWGR / pGR, 0., 1.);
        auto const xnWG = 1. - mass_mole_fractions.xnCG;
        ASSERT_NEAR(refence_xnWG, xnWG, 1.e-10);

        // The quotient of constituent partial densities and phase densities
        // must be equal to the mass fraction of those constituents in both
        // phases.
        ASSERT_NEAR(constituent_density.rho_W_GR / fluid_density.rho_GR,
                    1 - mass_mole_fractions.xmCG, 1.e-10);
        ASSERT_NEAR(rhoWLR() / fluid_density.rho_LR, mass_mole_fractions.xmWL,
                    1.e-10);

        ASSERT_NEAR(constituent_density.rho_C_GR / fluid_density.rho_GR,
                    mass_mole_fractions.xmCG, 1.e-10);
        ASSERT_NEAR(constituent_density.rho_C_LR / fluid_density.rho_LR,
                    1. - mass_mole_fractions.xmWL, 1.e-10);

        // Sum of constituent partial densities must be equal to phase
        // density
        ASSERT_NEAR(constituent_density.rho_C_GR + constituent_density.rho_W_GR,
                    fluid_density.rho_GR, 1.e-10);
        ASSERT_NEAR(constituent_density.rho_C_LR + rhoWLR(),
                    fluid_density.rho_LR, 1.e-10);

        // Liquid phase contains Water component only
        ASSERT_NEAR(1., mass_mole_fractions.xmWL, 1.e-10);

        // Gas density (ideal gas in this test):
        constexpr double R = MaterialLib::PhysicalConstant::IdealGasConstant;
        auto const MG =
            xnWG * molar_mass_water + mass_mole_fractions.xnCG * molar_mass_air;
        auto const rhoGR = pGR * MG / R / T;
        ASSERT_NEAR(rhoGR, fluid_density.rho_GR, 1.e-10);

        // Liquid density (linear EOS)
        auto const pLR = pGR - pCap;
        auto const rhoLR = rhoLR_ref * (1. + slope_pLR * (pLR - pLR_ref) +
                                        slope_T * (T - T_ref));
        ASSERT_NEAR(rhoLR, fluid_density.rho_LR, 1.e-10);

        // TODO: Test hAlpha, uAlpha
        // ASSERT_NEAR(hCG, cv.hCG, 1.0e-09);
        // ASSERT_NEAR(hWG, cv.hWG, 1.0e-09);
        // ASSERT_NEAR(hG, enthalpy.h_G, 1.0e-09);
        // ASSERT_NEAR(hL, enthalpy.h_L, 1.0e-09);
        // ASSERT_NEAR(uG, cv.uG, 1.0e-09);
        // ASSERT_NEAR(uL, cv.uL, 1.0e-09);

        /*auto const tortuosity = 1.;*/
        auto const diffusion_gas_phase = /*tortuosity * molar_mass_air *
                                         molar_mass_water / MG / MG * */
            diffusion_coefficient_vapour;
        ASSERT_NEAR(diffusion_gas_phase, cv.diffusion_coefficient_vapour,
                    1.0e-10);
    }
}

// Same test as above, but with constant gas density
TEST(ProcessLib, TH2MPhaseTransitionConstRho)
{
    namespace CR = ProcessLib::TH2M::ConstitutiveRelations;

    bool const density_is_constant = true;
    std::shared_ptr<MaterialPropertyLib::Medium> const& medium =
        Tests::createTestMaterial(MediumDefinition(density_is_constant));
    CR::MediaData media_data{*medium};

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media{
        {0, medium}};

    auto ptm = std::make_unique<CR::PhaseTransition>(media);

    auto const count = 200000;
    auto const pGR_min = 100000.;
    auto const pCap_min = 5000.;
    auto const T_min = 280.;

    auto const pGR_max = 200000.;
    auto const pCap_max = 5000000.;
    auto const T_max = 480.;

    auto const dpGR = (pGR_max - pGR_min) / count;
    auto const dpCap = (pCap_max - pCap_min) / count;
    auto const dT = (T_max - T_min) / count;
    auto const eps_pGR = pGR_max * 5.e-4;
    auto const eps_pCap = pCap_max * 2.e-5;
    auto const eps_T = T_max * 2.e-5;

    CR::PureLiquidDensityData rhoWLR;
    CR::PureLiquidDensityModel rhoWLR_model;

    CR::FluidEnthalpyData enthalpy;
    CR::MassMoleFractionsData mass_mole_fractions;
    CR::FluidDensityData fluid_density;
    CR::VapourPartialPressureData vapour_pressure;
    CR::ConstituentDensityData constituent_density;
    CR::PhaseTransitionData cv;
    CR::SpaceTimeData x_t{{},
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};

    for (std::size_t i = 0; i <= count; i++)
    {
        auto pGR = pGR_min + i * dpGR;
        auto pCap = pCap_min + i * dpCap;
        auto T = T_min + i * dT;
        CR::GasPressureData p_GR_data{pGR};
        CR::CapillaryPressureData p_cap_data{pCap};
        CR::TemperatureData T_data{T, T};

        // Perturb gas pressure
        rhoWLR_model.eval(x_t, media_data, CR::GasPressureData{pGR + eps_pGR},
                          p_cap_data, T_data, rhoWLR);
        ptm->eval(x_t, media_data, CR::GasPressureData{pGR + eps_pGR},
                  p_cap_data, T_data, rhoWLR, enthalpy, mass_mole_fractions,
                  fluid_density, vapour_pressure, constituent_density, cv);

        auto xmWG_plus = 1 - mass_mole_fractions.xmCG;
        auto rhoCGR_plus = constituent_density.rho_C_GR;
        auto rhoWGR_plus = constituent_density.rho_W_GR;

        rhoWLR_model.eval(x_t, media_data, CR::GasPressureData{pGR - eps_pGR},
                          p_cap_data, T_data, rhoWLR);
        ptm->eval(x_t, media_data, CR::GasPressureData{pGR - eps_pGR},
                  p_cap_data, T_data, rhoWLR, enthalpy, mass_mole_fractions,
                  fluid_density, vapour_pressure, constituent_density, cv);

        auto xmWG_minus = 1 - mass_mole_fractions.xmCG;
        auto rhoCGR_minus = constituent_density.rho_C_GR;
        auto rhoWGR_minus = constituent_density.rho_W_GR;

        // Unperturbed primary variables
        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        // Central difference derivatives
        auto const dxmWG_dpGR = (xmWG_plus - xmWG_minus) / (2. * eps_pGR);
        auto const drhoCGR_dpGR = (rhoCGR_plus - rhoCGR_minus) / (2. * eps_pGR);
        auto const drhoWGR_dpGR = (rhoWGR_plus - rhoWGR_minus) / (2. * eps_pGR);

        // Test composition derivatives w.r.t. gas pressure
        auto const tolerance_dpGR = 1e-4;

        if ((1. - xmWG_plus) > tolerance_dpGR)
        {
            ASSERT_NEAR(dxmWG_dpGR, cv.dxmWG_dpGR, tolerance_dpGR);
            ASSERT_NEAR(0., cv.drho_GR_dp_GR, tolerance_dpGR);
            ASSERT_NEAR(drhoCGR_dpGR, cv.drho_C_GR_dp_GR, tolerance_dpGR);
            ASSERT_NEAR(drhoWGR_dpGR, cv.drho_W_GR_dp_GR, tolerance_dpGR);
        }

        // Perturb capillary pressure
        rhoWLR_model.eval(x_t, media_data, p_GR_data,
                          CR::CapillaryPressureData{pCap + eps_pCap}, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data,
                  CR::CapillaryPressureData{pCap + eps_pCap}, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_plus = 1 - mass_mole_fractions.xmCG;
        rhoCGR_plus = constituent_density.rho_C_GR;
        rhoWGR_plus = constituent_density.rho_W_GR;

        rhoWLR_model.eval(x_t, media_data, p_GR_data,
                          CR::CapillaryPressureData{pCap - eps_pCap}, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data,
                  CR::CapillaryPressureData{pCap - eps_pCap}, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_minus = 1 - mass_mole_fractions.xmCG;
        rhoCGR_minus = constituent_density.rho_C_GR;
        rhoWGR_minus = constituent_density.rho_W_GR;

        // Unperturbed primary variables
        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        // Central difference derivatives
        auto const dxmWG_dpCap = (xmWG_plus - xmWG_minus) / (2. * eps_pCap);
        auto const drhoCGR_dpCap =
            (rhoCGR_plus - rhoCGR_minus) / (2. * eps_pCap);
        auto const drhoWGR_dpCap =
            (rhoWGR_plus - rhoWGR_minus) / (2. * eps_pCap);

        // Test composition derivatives w.r.t. capillary pressure
        auto const tolerance_dpCap = 1e-10;

        if ((1. - xmWG_plus) > tolerance_dpCap)
        {
            ASSERT_NEAR(dxmWG_dpCap, cv.dxmWG_dpCap, tolerance_dpCap);
            ASSERT_NEAR(0., cv.drho_GR_dp_cap, tolerance_dpCap);
            ASSERT_NEAR(drhoCGR_dpCap, cv.drho_C_GR_dp_cap, tolerance_dpCap);
            ASSERT_NEAR(drhoWGR_dpCap, cv.drho_W_GR_dp_cap, tolerance_dpCap);
        }

        // Perturb temperature
        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data,
                          CR::TemperatureData{T + eps_T, T}, rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data,
                  CR::TemperatureData{T + eps_T, T}, rhoWLR, enthalpy,
                  mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_plus = 1 - mass_mole_fractions.xmCG;
        rhoCGR_plus = constituent_density.rho_C_GR;
        rhoWGR_plus = constituent_density.rho_W_GR;

        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data,
                          CR::TemperatureData{T - eps_T, T}, rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data,
                  CR::TemperatureData{T - eps_T, T}, rhoWLR, enthalpy,
                  mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        xmWG_minus = 1 - mass_mole_fractions.xmCG;
        rhoCGR_minus = constituent_density.rho_C_GR;
        rhoWGR_minus = constituent_density.rho_W_GR;

        // Unperturbed primary variables
        rhoWLR_model.eval(x_t, media_data, p_GR_data, p_cap_data, T_data,
                          rhoWLR);
        ptm->eval(x_t, media_data, p_GR_data, p_cap_data, T_data, rhoWLR,
                  enthalpy, mass_mole_fractions, fluid_density, vapour_pressure,
                  constituent_density, cv);

        // Central difference derivatives
        auto const dxmWG_dT = (xmWG_plus - xmWG_minus) / (2. * eps_T);
        auto const drhoCGR_dT = (rhoCGR_plus - rhoCGR_minus) / (2. * eps_T);
        auto const drhoWGR_dT = (rhoWGR_plus - rhoWGR_minus) / (2. * eps_T);

        // Test composition derivatives w.r.t. temperature
        auto const tolerance_dT = 1.e-4;

        if ((1. - xmWG_plus) > tolerance_dT)
        {
            ASSERT_NEAR(dxmWG_dT, cv.dxmWG_dT, tolerance_dT);
            ASSERT_NEAR(0., cv.drho_GR_dT, tolerance_dT);
            ASSERT_NEAR(drhoCGR_dT, cv.drho_C_GR_dT, tolerance_dT);
            ASSERT_NEAR(drhoWGR_dT, cv.drho_W_GR_dT, tolerance_dT);
        }
        // Test mixture composition
        // Vapour mole fraction is the quotient of vapour pressure and gas
        // phase pressure
        auto const refence_xnWG =
            std::clamp(vapour_pressure.pWGR / pGR, 0., 1.);
        auto const xnWG = 1. - mass_mole_fractions.xnCG;
        ASSERT_NEAR(refence_xnWG, xnWG, 1.e-10);

        // The quotient of constituent partial densities and phase densities
        // must be equal to the mass fraction of those constituents in both
        // phases.
        ASSERT_NEAR(constituent_density.rho_W_GR / fluid_density.rho_GR,
                    1 - mass_mole_fractions.xmCG, 1.e-10);
        ASSERT_NEAR(rhoWLR() / fluid_density.rho_LR, mass_mole_fractions.xmWL,
                    1.e-10);

        ASSERT_NEAR(constituent_density.rho_C_GR / fluid_density.rho_GR,
                    mass_mole_fractions.xmCG, 1.e-10);
        ASSERT_NEAR(constituent_density.rho_C_LR / fluid_density.rho_LR,
                    1. - mass_mole_fractions.xmWL, 1.e-10);

        // Sum of constituent partial densities must be equal to phase
        // density
        ASSERT_NEAR(constituent_density.rho_C_GR + constituent_density.rho_W_GR,
                    fluid_density.rho_GR, 1.e-10);
        ASSERT_NEAR(constituent_density.rho_C_LR + rhoWLR(),
                    fluid_density.rho_LR, 1.e-10);

        // Liquid phase contains Water component only
        ASSERT_NEAR(1., mass_mole_fractions.xmWL, 1.e-10);

        // Gas density
        auto const rhoGR = constant_gas_density;
        ASSERT_NEAR(rhoGR, fluid_density.rho_GR, 1.e-10);

        // Liquid density (linear EOS)
        auto const pLR = pGR - pCap;
        auto const rhoLR = rhoLR_ref * (1. + slope_pLR * (pLR - pLR_ref) +
                                        slope_T * (T - T_ref));
        ASSERT_NEAR(rhoLR, fluid_density.rho_LR, 1.e-10);

        // TODO: Test hAlpha, uAlpha
        // ASSERT_NEAR(hCG, cv.hCG, 1.0e-09);
        // ASSERT_NEAR(hWG, cv.hWG, 1.0e-09);
        // ASSERT_NEAR(hG, enthalpy.h_G, 1.0e-09);
        // ASSERT_NEAR(hL, enthalpy.h_L, 1.0e-09);
        // ASSERT_NEAR(uG, cv.uG, 1.0e-09);
        // ASSERT_NEAR(uL, cv.uL, 1.0e-09);

        /*auto const tortuosity = 1.;*/
        auto const diffusion_gas_phase = /*tortuosity * molar_mass_air *
                                         molar_mass_water / MG / MG * */
            diffusion_coefficient_vapour;
        ASSERT_NEAR(diffusion_gas_phase, cv.diffusion_coefficient_vapour,
                    1.0e-10);
    }
}
