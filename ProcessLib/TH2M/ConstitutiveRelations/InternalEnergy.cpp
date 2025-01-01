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
void InternalEnergyModel::eval(FluidDensityData const& fluid_density_data,
                               PhaseTransitionData const& phase_transition_data,
                               PorosityData const& porosity_data,
                               SaturationData const& S_L_data,
                               SolidDensityData const& solid_density_data,
                               SolidEnthalpyData const& solid_enthalpy_data,
                               InternalEnergyData& internal_energy_data) const
{
    auto const phi_L = S_L_data.S_L * porosity_data.phi;
    auto const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;
    double const phi_S = 1. - porosity_data.phi;

    auto const u_S = solid_enthalpy_data.h_S;

    internal_energy_data() =
        phi_G * fluid_density_data.rho_GR * phase_transition_data.uG +
        phi_L * fluid_density_data.rho_LR * phase_transition_data.uL +
        phi_S * solid_density_data.rho_SR * u_S;
}

void InternalEnergyModel::dEval(
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    PorosityDerivativeData const& porosity_d_data,
    SaturationData const& S_L_data,
    SolidDensityData const& solid_density_data,
    SolidDensityDerivativeData const& solid_density_d_data,
    SolidEnthalpyData const& solid_enthalpy_data,
    SolidHeatCapacityData const& solid_heat_capacity_data,
    EffectiveVolumetricInternalEnergyDerivatives&
        effective_volumetric_internal_energy_d_data) const
{
    auto const phi_L = S_L_data.S_L * porosity_data.phi;
    auto const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;
    double const phi_S = 1. - porosity_data.phi;

    auto const u_S = solid_enthalpy_data.h_S;
    effective_volumetric_internal_energy_d_data.drho_u_eff_dT =
        phi_G * phase_transition_data.drho_GR_dT * phase_transition_data.uG +
        phi_G * fluid_density_data.rho_GR * phase_transition_data.du_G_dT +
        phi_L * phase_transition_data.drho_LR_dT * phase_transition_data.uL +
        phi_L * fluid_density_data.rho_LR * phase_transition_data.du_L_dT +
        phi_S * solid_density_d_data.drho_SR_dT * u_S +
        phi_S * solid_density_data.rho_SR * solid_heat_capacity_data() -
        porosity_d_data.dphi_dT * solid_density_data.rho_SR * u_S;

    // From p_LR = p_GR - p_cap it follows for
    // drho_LR/dp_GR = drho_LR/dp_LR * dp_LR/dp_GR
    //               = drho_LR/dp_LR * (dp_GR/dp_GR - dp_cap/dp_GR)
    //               = drho_LR/dp_LR * (1 - 0)
    double const drho_LR_dp_GR = phase_transition_data.drho_LR_dp_LR;
    double const drho_LR_dp_cap = -phase_transition_data.drho_LR_dp_LR;
    // drho_GR_dp_cap = 0;

    effective_volumetric_internal_energy_d_data.drho_u_eff_dp_GR =
        /*(dphi_G_dp_GR = 0) * fluid_density_data.rho_GR *
           phase_transition_data.uG +*/
        phi_G * phase_transition_data.drho_GR_dp_GR * phase_transition_data.uG +
        phi_G * fluid_density_data.rho_GR * phase_transition_data.du_G_dp_GR +
        /*(dphi_L_dp_GR = 0) * fluid_density_data.rho_LR *
           phase_transition_data.uL +*/
        phi_L * drho_LR_dp_GR * phase_transition_data.uL +
        phi_L * fluid_density_data.rho_LR * phase_transition_data.du_L_dp_GR;

    effective_volumetric_internal_energy_d_data.drho_u_eff_dp_cap =
        -porosity_d_data.dphi_L_dp_cap * fluid_density_data.rho_GR *
            phase_transition_data.uG +
        /*phi_G * (drho_GR_dp_cap = 0) * phase_transition_data.uG +*/
        porosity_d_data.dphi_L_dp_cap * fluid_density_data.rho_LR *
            phase_transition_data.uL +
        phi_L * drho_LR_dp_cap * phase_transition_data.uL +
        phi_L * fluid_density_data.rho_LR * phase_transition_data.du_L_dp_cap;
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
