/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "InternalEnergy.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
void EffectiveVolumetricEnthalpyModel::eval(
    FluidDensityData const& fluid_density_data,
    FluidEnthalpyData const& fluid_enthalpy_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    SolidDensityData const& solid_density_data,
    SolidEnthalpyData const& solid_enthalpy_data,
    EffectiveVolumetricEnthalpy& effective_volumetric_enthalpy_data) const
{
    auto const phi_L = S_L_data.S_L * porosity_data.phi;
    auto const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;
    double const phi_S = 1. - porosity_data.phi;

    effective_volumetric_enthalpy_data.rho_h_eff =
        phi_G * fluid_density_data.rho_GR * fluid_enthalpy_data.h_G +
        phi_L * fluid_density_data.rho_LR * fluid_enthalpy_data.h_L +
        phi_S * solid_density_data.rho_SR * solid_enthalpy_data.h_S;
}

void EffectiveVolumetricEnthalpyModel::dEval(
    FluidDensityData const& fluid_density_data,
    FluidEnthalpyData const& fluid_enthalpy_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    PorosityDerivativeData const& porosity_d_data,
    SaturationData const& S_L_data,
    SolidDensityData const& solid_density_data,
    SolidDensityDerivativeData const& solid_density_d_data,
    SolidEnthalpyData const& solid_enthalpy_data,
    SolidHeatCapacityData const& solid_heat_capacity_data,
    EffectiveVolumetricEnthalpyDerivatives&
        effective_volumetric_enthalpy_d_data) const
{
    auto const phi_L = S_L_data.S_L * porosity_data.phi;
    auto const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;
    double const phi_S = 1. - porosity_data.phi;

    // From p_LR = p_GR - p_cap it follows for
    // drho_LR/dp_GR = drho_LR/dp_LR * dp_LR/dp_GR
    //               = drho_LR/dp_LR * (dp_GR/dp_GR - dp_cap/dp_GR)
    //               = drho_LR/dp_LR * (1 - 0)
    double const drho_LR_dp_GR = phase_transition_data.drho_LR_dp_LR;
    double const drho_LR_dp_cap = -phase_transition_data.drho_LR_dp_LR;
    // drho_GR_dp_cap = 0;

    effective_volumetric_enthalpy_d_data.drho_h_eff_dp_GR =
        /*(dphi_G_dp_GR = 0) * fluid_density_data.rho_GR *
            fluid_enthalpy_data.h_G +*/
        phi_G * phase_transition_data.drho_GR_dp_GR * fluid_enthalpy_data.h_G +
        /*(dphi_L_dp_GR = 0) * fluid_density_data.rho_LR *
            fluid_enthalpy_data.h_L +*/
        phi_L * drho_LR_dp_GR * fluid_enthalpy_data.h_L;
    effective_volumetric_enthalpy_d_data.drho_h_eff_dp_cap =
        porosity_d_data.dphi_L_dp_cap * fluid_density_data.rho_GR *
            fluid_enthalpy_data.h_G +
        /*phi_G * (drho_GR_dp_cap = 0) * fluid_enthalpy_data.h_G +*/
        porosity_d_data.dphi_L_dp_cap * fluid_density_data.rho_LR *
            fluid_enthalpy_data.h_L +
        phi_L * drho_LR_dp_cap * fluid_enthalpy_data.h_L;

    // TODO (naumov) Extend for temperature dependent porosities.
    constexpr double dphi_G_dT = 0;
    constexpr double dphi_L_dT = 0;
    effective_volumetric_enthalpy_d_data.drho_h_eff_dT =
        dphi_G_dT * fluid_density_data.rho_GR * fluid_enthalpy_data.h_G +
        phi_G * phase_transition_data.drho_GR_dT * fluid_enthalpy_data.h_G +
        phi_G * fluid_density_data.rho_GR * phase_transition_data.dh_G_dT +
        dphi_L_dT * fluid_density_data.rho_LR * fluid_enthalpy_data.h_L +
        phi_L * phase_transition_data.drho_LR_dT * fluid_enthalpy_data.h_L +
        phi_L * fluid_density_data.rho_LR * phase_transition_data.dh_L_dT -
        porosity_d_data.dphi_dT * solid_density_data.rho_SR *
            solid_enthalpy_data.h_S +
        phi_S * solid_density_d_data.drho_SR_dT * solid_enthalpy_data.h_S +
        phi_S * solid_density_data.rho_SR * solid_heat_capacity_data();
}

void SolidEnthalpyModel::eval(
    SolidHeatCapacityData const& solid_heat_capacity_data,
    TemperatureData const& T_data,
    SolidEnthalpyData& solid_enthalpy_data) const
{
    solid_enthalpy_data.h_S = solid_heat_capacity_data() * T_data.T;
}

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
