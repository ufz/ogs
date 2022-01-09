/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseTransitionNone.h"

#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
namespace TH2M
{
PhaseTransitionNone::PhaseTransitionNone(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
    : PhaseTransitionModel(media)
{
    DBUG("Create PhaseTransitionNone constitutive model.");

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

PhaseTransitionModelVariables PhaseTransitionNone::updateConstitutiveVariables(
    PhaseTransitionModelVariables const& phase_transition_model_variables,
    const MaterialPropertyLib::Medium* medium,
    MaterialPropertyLib::VariableArray variables,
    ParameterLib::SpatialPosition pos, double const t, double const dt) const
{
    // primary variables
    auto const pGR = std::get<double>(variables[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)]);
    auto const T = std::get<double>(variables[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)]);

    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& gas_phase = medium->phase("Gas");

    // copy previous state before modification.
    PhaseTransitionModelVariables cv = phase_transition_model_variables;

    // C-component is only component in the gas phase
    cv.xnCG = 1.;
    cv.xmCG = 1.;

    auto const M =
        gas_phase.property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, pos, t, dt);

    variables[static_cast<int>(MaterialPropertyLib::Variable::molar_mass)] = M;

    cv.rhoGR = gas_phase.property(MaterialPropertyLib::PropertyType::density)
                   .template value<double>(variables, pos, t, dt);
    cv.muGR = gas_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                  .template value<double>(variables, pos, t, dt);
    cv.lambdaGR =
        gas_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variables, pos, t, dt);

    cv.rhoCGR = cv.rhoGR;

    // W-component is only component in the liquid phase
    cv.xmWL = 1.;

    cv.rhoLR = liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                   .template value<double>(variables, pos, t, dt);

    cv.muLR =
        liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
            .template value<double>(variables, pos, t, dt);

    cv.lambdaLR =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variables, pos, t, dt);

    cv.rhoWLR = cv.rhoLR;

    // specific heat capacities
    auto const cpG =
        gas_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);

    auto const cpL =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, pos, t, dt);

    // specific phase enthalpies
    cv.hG = cpG * T;
    cv.hL = cpL * T;
    cv.dh_G_dT = cpG;
    cv.dh_L_dT = cpL;

    // specific inner energies
    cv.uG = cv.hG - pGR / cv.rhoGR;
    cv.uL = cv.hL;

    auto const drho_GR_dT =
        gas_phase[MaterialPropertyLib::PropertyType::density]
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, t, dt);
    cv.du_G_dT = cpG + pGR * drho_GR_dT / cv.rhoGR / cv.rhoGR;

    cv.du_L_dT = cpL;

    cv.drho_GR_dp_GR =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::phase_pressure, pos,
                t, dt);
    cv.drho_LR_dp_LR =
        liquid_phase[MaterialPropertyLib::PropertyType::density]
            .template dValue<double>(
                variables, MaterialPropertyLib::Variable::liquid_phase_pressure,
                pos, t, dt);
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
                                     pos, t, dt);
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
    return cv;
}
}  // namespace TH2M
}  // namespace ProcessLib
