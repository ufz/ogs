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
#include "MaterialLib/PhysicalConstant.h"
#include "ProcessLib/TH2M/PhaseTransitionModels/PhaseTransitionEvaporation.h"
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
static const double viscosity_air = 1.e-5;
static const double viscosity_water = 1.e-3;

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
    m << Tests::makeConstantPropertyElement("viscosity", viscosity_air);
    m << Tests::makeConstantPropertyElement("thermal_conductivity",
                                            thermal_conductivity_air);

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

    // liquid phase properties
    m << "<properties>\n";
    m << Tests::makeConstantPropertyElement("viscosity", viscosity_water);
    m << Tests::makeConstantPropertyElement("specific_heat_capacity",
                                            specific_heat_capacity_water);

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
    m << "</phases> </medium>\n";

    return m.str();
}

TEST(ProcessLib, TH2MPhaseTransitionEvaporation)
{
    bool const density_is_constant = false;
    std::shared_ptr<MaterialPropertyLib::Medium> const& medium =
        Tests::createTestMaterial(MediumDefinition(density_is_constant));

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media{
        {0, medium}};

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    auto ptm =
        std::make_unique<ProcessLib::TH2M::PhaseTransitionEvaporation>(media);

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

    for (std::size_t i = 0; i <= count; i++)
    {
        auto pGR = pGR_min + i * dpGR;
        auto pCap = pCap_min + i * dpCap;
        auto T = T_min + i * dT;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        auto const& cv = ptm->cv;

        // Perturb gas pressure
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR + eps_pGR;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        auto xmWG_plus = ptm->cv.xmWG;
        auto rhoGR_plus = ptm->cv.rhoGR;
        auto rhoCGR_plus = ptm->cv.rhoCGR;
        auto rhoWGR_plus = ptm->cv.rhoWGR;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR - eps_pGR;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        auto xmWG_minus = ptm->cv.xmWG;
        auto rhoGR_minus = ptm->cv.rhoGR;
        auto rhoCGR_minus = ptm->cv.rhoCGR;
        auto rhoWGR_minus = ptm->cv.rhoWGR;

        // Unperturbed primary variables
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

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
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] =
            pCap + eps_pCap;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_plus = ptm->cv.xmWG;
        rhoGR_plus = ptm->cv.rhoGR;
        rhoCGR_plus = ptm->cv.rhoCGR;
        rhoWGR_plus = ptm->cv.rhoWGR;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] =
            pCap - eps_pCap;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_minus = ptm->cv.xmWG;
        rhoGR_minus = ptm->cv.rhoGR;
        rhoCGR_minus = ptm->cv.rhoCGR;
        rhoWGR_minus = ptm->cv.rhoWGR;

        // Unperturbed primary variables
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

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
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T + eps_T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_plus = ptm->cv.xmWG;
        rhoGR_plus = ptm->cv.rhoGR;
        rhoCGR_plus = ptm->cv.rhoCGR;
        rhoWGR_plus = ptm->cv.rhoWGR;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T - eps_T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_minus = ptm->cv.xmWG;
        rhoGR_minus = ptm->cv.rhoGR;
        rhoCGR_minus = ptm->cv.rhoCGR;
        rhoWGR_minus = ptm->cv.rhoWGR;

        // Unperturbed primary variables
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

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
        auto const refence_xnWG = std::clamp(cv.pWGR / pGR, 0., 1.);
        ASSERT_NEAR(refence_xnWG, cv.xnWG, 1.e-10);

        // The quotient of constituent partial densities and phase densities
        // must be equal to the mass fraction of those constituents in both
        // phases.
        ASSERT_NEAR(cv.rhoWGR / cv.rhoGR, cv.xmWG, 1.e-10);
        ASSERT_NEAR(cv.rhoWLR / cv.rhoLR, cv.xmWL, 1.e-10);

        ASSERT_NEAR(cv.rhoCGR / cv.rhoGR, 1. - cv.xmWG, 1.e-10);
        ASSERT_NEAR(cv.rhoCLR / cv.rhoLR, 1. - cv.xmWL, 1.e-10);

        // Sum of constituent partial densities must be equal to phase
        // density
        ASSERT_NEAR(cv.rhoCGR + cv.rhoWGR, cv.rhoGR, 1.e-10);
        ASSERT_NEAR(cv.rhoCLR + cv.rhoWLR, cv.rhoLR, 1.e-10);

        // Liquid phase contains Water component only
        ASSERT_NEAR(1., cv.xmWL, 1.e-10);

        // Gas density (ideal gas in this test):
        constexpr double R = MaterialLib::PhysicalConstant::IdealGasConstant;
        auto const xnCG = 1. - cv.xnWG;
        auto const MG = cv.xnWG * molar_mass_water + xnCG * molar_mass_air;
        auto const rhoGR = pGR * MG / R / T;
        ASSERT_NEAR(rhoGR, cv.rhoGR, 1.e-10);

        // Liquid density (linear EOS)
        auto const pLR = pGR - pCap;
        auto const rhoLR = rhoLR_ref * (1. + slope_pLR * (pLR - pLR_ref) +
                                        slope_T * (T - T_ref));
        ASSERT_NEAR(rhoLR, cv.rhoLR, 1.e-10);

        // TODO: Test hAlpha, uAlpha
        // ASSERT_NEAR(hCG, cv.hCG, 1.0e-09);
        // ASSERT_NEAR(hWG, cv.hWG, 1.0e-09);
        // ASSERT_NEAR(hG, cv.hG, 1.0e-09);
        // ASSERT_NEAR(hL, cv.hL, 1.0e-09);
        // ASSERT_NEAR(uG, cv.uG, 1.0e-09);
        // ASSERT_NEAR(uL, cv.uL, 1.0e-09);

        /*auto const tortuosity = 1.;*/
        auto const diffusion_gas_phase = /*tortuosity * molar_mass_air *
                                         molar_mass_water / MG / MG * */
            diffusion_coefficient_vapour;
        ASSERT_NEAR(diffusion_gas_phase, cv.diffusion_coefficient_vapour,
                    1.0e-10);

        ASSERT_NEAR(viscosity_air, cv.muGR, 1.0e-10);
        ASSERT_NEAR(viscosity_water, cv.muLR, 1.0e-10);
    }
}

// Same test as above, but with constant gas density
TEST(ProcessLib, TH2MPhaseTransitionEvaporationConstRho)
{
    bool const density_is_constant = true;
    std::shared_ptr<MaterialPropertyLib::Medium> const& medium =
        Tests::createTestMaterial(MediumDefinition(density_is_constant));

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media{
        {0, medium}};

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    auto ptm =
        std::make_unique<ProcessLib::TH2M::PhaseTransitionEvaporation>(media);

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

    for (std::size_t i = 0; i <= count; i++)
    {
        auto pGR = pGR_min + i * dpGR;
        auto pCap = pCap_min + i * dpCap;
        auto T = T_min + i * dT;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        auto const& cv = ptm->cv;

        // Perturb gas pressure
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR + eps_pGR;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        auto xmWG_plus = ptm->cv.xmWG;
        auto rhoCGR_plus = ptm->cv.rhoCGR;
        auto rhoWGR_plus = ptm->cv.rhoWGR;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR - eps_pGR;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        auto xmWG_minus = ptm->cv.xmWG;
        auto rhoCGR_minus = ptm->cv.rhoCGR;
        auto rhoWGR_minus = ptm->cv.rhoWGR;

        // Unperturbed primary variables
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

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
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] =
            pCap + eps_pCap;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_plus = ptm->cv.xmWG;
        rhoCGR_plus = ptm->cv.rhoCGR;
        rhoWGR_plus = ptm->cv.rhoWGR;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] =
            pCap - eps_pCap;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_minus = ptm->cv.xmWG;
        rhoCGR_minus = ptm->cv.rhoCGR;
        rhoWGR_minus = ptm->cv.rhoWGR;

        // Unperturbed primary variables
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

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
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T + eps_T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_plus = ptm->cv.xmWG;
        rhoCGR_plus = ptm->cv.rhoCGR;
        rhoWGR_plus = ptm->cv.rhoWGR;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T - eps_T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

        xmWG_minus = ptm->cv.xmWG;
        rhoCGR_minus = ptm->cv.rhoCGR;
        rhoWGR_minus = ptm->cv.rhoWGR;

        // Unperturbed primary variables
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = pGR;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pCap;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        ptm->computeConstitutiveVariables(medium.get(), variable_array, pos,
                                          time, dt);

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
        auto const refence_xnWG = std::clamp(cv.pWGR / pGR, 0., 1.);
        ASSERT_NEAR(refence_xnWG, cv.xnWG, 1.e-10);

        // The quotient of constituent partial densities and phase densities
        // must be equal to the mass fraction of those constituents in both
        // phases.
        ASSERT_NEAR(cv.rhoWGR / cv.rhoGR, cv.xmWG, 1.e-10);
        ASSERT_NEAR(cv.rhoWLR / cv.rhoLR, cv.xmWL, 1.e-10);

        ASSERT_NEAR(cv.rhoCGR / cv.rhoGR, 1. - cv.xmWG, 1.e-10);
        ASSERT_NEAR(cv.rhoCLR / cv.rhoLR, 1. - cv.xmWL, 1.e-10);

        // Sum of constituent partial densities must be equal to phase
        // density
        ASSERT_NEAR(cv.rhoCGR + cv.rhoWGR, cv.rhoGR, 1.e-10);
        ASSERT_NEAR(cv.rhoCLR + cv.rhoWLR, cv.rhoLR, 1.e-10);

        // Liquid phase contains Water component only
        ASSERT_NEAR(1., cv.xmWL, 1.e-10);

        // Gas density
        auto const rhoGR = constant_gas_density;
        ASSERT_NEAR(rhoGR, cv.rhoGR, 1.e-10);

        // Liquid density (linear EOS)
        auto const pLR = pGR - pCap;
        auto const rhoLR = rhoLR_ref * (1. + slope_pLR * (pLR - pLR_ref) +
                                        slope_T * (T - T_ref));
        ASSERT_NEAR(rhoLR, cv.rhoLR, 1.e-10);

        // TODO: Test hAlpha, uAlpha
        // ASSERT_NEAR(hCG, cv.hCG, 1.0e-09);
        // ASSERT_NEAR(hWG, cv.hWG, 1.0e-09);
        // ASSERT_NEAR(hG, cv.hG, 1.0e-09);
        // ASSERT_NEAR(hL, cv.hL, 1.0e-09);
        // ASSERT_NEAR(uG, cv.uG, 1.0e-09);
        // ASSERT_NEAR(uL, cv.uL, 1.0e-09);

        /*auto const tortuosity = 1.;*/
        auto const diffusion_gas_phase = /*tortuosity * molar_mass_air *
                                         molar_mass_water / MG / MG * */
            diffusion_coefficient_vapour;
        ASSERT_NEAR(diffusion_gas_phase, cv.diffusion_coefficient_vapour,
                    1.0e-10);

        ASSERT_NEAR(viscosity_air, cv.muGR, 1.0e-10);
        ASSERT_NEAR(viscosity_water, cv.muLR, 1.0e-10);
    }
}