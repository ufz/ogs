/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "NoPhaseTransition.h"

#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
NoPhaseTransition::NoPhaseTransition(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
    : PhaseTransitionModel(media)
{
    DBUG("Create NoPhaseTransition constitutive model.");

    // check for minimum requirement definitions in media object
    std::array const required_gas_properties = {
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::molar_mass};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::specific_heat_capacity};

    for (auto const& m : media)
    {
        checkRequiredProperties(m.second->phase("Gas"),
                                required_gas_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
    }
}

void NoPhaseTransition::eval(SpaceTimeData const& x_t,
                             MediaData const& media_data,
                             GasPressureData const& p_GR,
                             CapillaryPressureData const& p_cap,
                             TemperatureData const& T_data,
                             PureLiquidDensityData const& rho_W_LR,
                             ViscosityData& viscosity_data,
                             EnthalpyData& enthalpy_data,
                             MassMoleFractionsData& mass_mole_fractions_data,
                             PhaseTransitionData& cv) const
{
    MaterialPropertyLib::VariableArray variables;

    // primary variables
    auto const pGR = p_GR();
    auto const pCap = p_cap();
    auto const T = T_data.T;
    variables.gas_phase_pressure = pGR;
    variables.temperature = T;

    auto const& liquid_phase = media_data.liquid_phase;
    auto const& gas_phase = media_data.gas_phase;

    // C-component is only component in the gas phase
    cv.xnWG = 0.;
    cv.xmWG = 0.;
    mass_mole_fractions_data.xnCG = 1. - cv.xnWG;
    mass_mole_fractions_data.xmCG = 1. - cv.xmWG;

    auto const M =
        gas_phase.property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    variables.molar_mass = M;

    cv.rhoGR = gas_phase.property(MaterialPropertyLib::PropertyType::density)
                   .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    viscosity_data.mu_GR =
        gas_phase.property(MaterialPropertyLib::PropertyType::viscosity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    cv.rhoCGR = cv.rhoGR;

    // W-component is only component in the liquid phase
    mass_mole_fractions_data.xmWL = 1.;

    auto const pLR = pGR - pCap;
    variables.liquid_phase_pressure = pLR;
    cv.rhoLR = rho_W_LR();
    variables.density = cv.rhoLR;

    viscosity_data.mu_LR =
        liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // specific heat capacities
    auto const cpG =
        gas_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    auto const cpL =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    // specific phase enthalpies
    enthalpy_data.h_G = cpG * T;
    enthalpy_data.h_L = cpL * T;
    cv.dh_G_dT = cpG;
    cv.dh_L_dT = cpL;

    // specific inner energies
    cv.uG = enthalpy_data.h_G - pGR / cv.rhoGR;
    cv.uL = enthalpy_data.h_L;

    auto const drho_GR_dT =
        gas_phase[MaterialPropertyLib::PropertyType::density]
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);
    cv.du_G_dT = cpG + pGR * drho_GR_dT / cv.rhoGR / cv.rhoGR;

    cv.du_L_dT = cpL;

    cv.drho_GR_dp_GR =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::gas_phase_pressure,
                x_t.x, x_t.t, x_t.dt);
    cv.drho_LR_dp_LR =
        liquid_phase[MaterialPropertyLib::PropertyType::density]
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::liquid_phase_pressure,
                x_t.x, x_t.t, x_t.dt);
    cv.drho_LR_dp_GR = cv.drho_LR_dp_LR;

    cv.du_G_dp_GR =
        -1 / cv.rhoGR + pGR * cv.drho_GR_dp_GR / cv.rhoGR / cv.rhoGR;

    cv.drho_C_GR_dp_GR = cv.drho_GR_dp_GR;
    cv.drho_C_LR_dp_LR = 0;
    cv.drho_C_LR_dp_GR = 0;
    cv.drho_C_GR_dT = drho_GR_dT;

    auto const drho_LR_dT =
        liquid_phase[MaterialPropertyLib::PropertyType::density]
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);
    cv.drho_C_LR_dT = 0;

    cv.du_L_dp_GR = 0;
    cv.du_L_dp_cap = 0;
    /* TODO update to the following when uL has same structure as the uG:
        +-1 / cv.rhoLR + pLR * cv.drho_LR_dp_cap / cv.rhoLR / cv.rhoLR;
    */

    cv.drho_W_LR_dp_LR = cv.drho_LR_dp_LR;
    cv.drho_W_LR_dp_GR = cv.drho_LR_dp_GR;
    cv.drho_W_LR_dT = drho_LR_dT;
    cv.drho_W_GR_dT = 0;
    cv.drho_W_GR_dp_GR = 0;
    cv.drho_W_GR_dp_cap = 0;
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
