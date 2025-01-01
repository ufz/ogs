/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
                             FluidEnthalpyData& fluid_enthalpy_data,
                             MassMoleFractionsData& mass_mole_fractions_data,
                             FluidDensityData& fluid_density_data,
                             VapourPartialPressureData& vapour_pressure_data,
                             ConstituentDensityData& constituent_density_data,
                             PhaseTransitionData& cv) const
{
    MaterialPropertyLib::VariableArray variables;

    // primary variables
    auto const pGR = p_GR();
    auto const pCap = p_cap();
    auto const T = T_data.T;
    variables.gas_phase_pressure = pGR;
    variables.temperature = T;

    auto const& liquid_phase = media_data.liquid;
    auto const& gas_phase = media_data.gas;

    vapour_pressure_data.pWGR = 0;

    // C-component is only component in the gas phase
    cv.dxmWG_dpGR = 0.;
    cv.dxmWG_dpCap = 0.;
    cv.dxmWG_dT = 0.;
    mass_mole_fractions_data.xnCG = 1.;
    mass_mole_fractions_data.xmCG = 1.;

    auto const M =
        gas_phase.property(MaterialPropertyLib::PropertyType::molar_mass)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    variables.molar_mass = M;

    fluid_density_data.rho_GR =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    constituent_density_data.rho_C_GR = fluid_density_data.rho_GR;
    constituent_density_data.rho_W_GR = 0;
    constituent_density_data.rho_C_LR = 0;

    // W-component is only component in the liquid phase
    mass_mole_fractions_data.xmWL = 1.;
    cv.dxmWL_dpGR = 0;
    cv.dxmWL_dpCap = 0;
    cv.dxmWL_dT = 0;

    auto const pLR = pGR - pCap;
    variables.liquid_phase_pressure = pLR;
    fluid_density_data.rho_LR = rho_W_LR();
    variables.density = fluid_density_data.rho_LR;

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
    fluid_enthalpy_data.h_G = cpG * T;
    fluid_enthalpy_data.h_L = cpL * T;
    cv.dh_G_dT = cpG;
    cv.dh_L_dT = cpL;

    // specific inner energies
    cv.uG = fluid_enthalpy_data.h_G - pGR / fluid_density_data.rho_GR;
    cv.uL = fluid_enthalpy_data.h_L;

    cv.drho_GR_dT =
        gas_phase[MaterialPropertyLib::PropertyType::density]
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);
    cv.du_G_dT = cpG + pGR * cv.drho_GR_dT / fluid_density_data.rho_GR /
                           fluid_density_data.rho_GR;

    cv.du_L_dT = cpL;

    cv.drho_GR_dp_cap = 0;

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

    cv.du_G_dp_GR = -1 / fluid_density_data.rho_GR +
                    pGR * cv.drho_GR_dp_GR / fluid_density_data.rho_GR /
                        fluid_density_data.rho_GR;

    cv.drho_C_GR_dp_GR = cv.drho_GR_dp_GR;
    cv.drho_C_GR_dp_cap = 0;
    cv.drho_C_GR_dT = cv.drho_GR_dT;

    cv.drho_LR_dT =
        liquid_phase[MaterialPropertyLib::PropertyType::density]
            .template dValue<double>(variables,
                                     MaterialPropertyLib::Variable::temperature,
                                     x_t.x, x_t.t, x_t.dt);

    cv.drho_W_LR_dp_LR = cv.drho_LR_dp_LR;
    cv.drho_W_LR_dp_GR = cv.drho_LR_dp_LR;
    cv.drho_W_LR_dT = cv.drho_LR_dT;
    cv.drho_W_GR_dT = 0;
    cv.drho_W_GR_dp_GR = 0;
    cv.drho_W_GR_dp_cap = 0;

    cv.diffusion_coefficient_vapour = 0.;
    cv.diffusion_coefficient_solute = 0.;

    cv.hCG = 0;
    cv.hWG = 0;
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
