/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
void FT1Model::eval(
    double const dt,
    InternalEnergyData const& internal_energy_data,
    PrevState<InternalEnergyData> const& internal_energy_data_prev,
    FT1Data& fT_1) const
{
    if (dt == 0)
    {
        fT_1.m = 0;
        return;
    }

    auto const rho_u_eff_dot =
        (internal_energy_data() - **internal_energy_data_prev) / dt;
    fT_1.m = rho_u_eff_dot;
}

void FT1Model::dEval(double const dt,
                     EffectiveVolumetricInternalEnergyDerivatives const&
                         effective_volumetric_internal_energy_d_data,
                     FT1DerivativeData& dfT_1) const
{
    if (dt == 0)
    {
        dfT_1.dp_GR = 0;
        dfT_1.dp_cap = 0;
        dfT_1.dT = 0;
        return;
    }

    dfT_1.dp_GR =
        effective_volumetric_internal_energy_d_data.drho_u_eff_dp_GR / dt;

    dfT_1.dp_cap =
        effective_volumetric_internal_energy_d_data.drho_u_eff_dp_cap / dt;

    dfT_1.dT = effective_volumetric_internal_energy_d_data.drho_u_eff_dT / dt;
}

template <int DisplacementDim>
void FT2Model<DisplacementDim>::eval(
    DarcyVelocityData<DisplacementDim> const& darcy_velocity_data,
    FluidDensityData const& fluid_density_data,
    FluidEnthalpyData const& fluid_enthalpy_data,
    FT2Data<DisplacementDim>& fT_2) const
{
    fT_2.A.noalias() = fluid_density_data.rho_GR * fluid_enthalpy_data.h_G *
                           darcy_velocity_data.w_GS +
                       fluid_density_data.rho_LR * fluid_enthalpy_data.h_L *
                           darcy_velocity_data.w_LS;
}

template <int DisplacementDim>
void FT2Model<DisplacementDim>::dEval(
    DarcyVelocityData<DisplacementDim> const& darcy_velocity_data,
    FluidDensityData const& fluid_density_data,
    FluidEnthalpyData const& fluid_enthalpy_data,
    PermeabilityData<DisplacementDim> const& permeability_data,
    PhaseTransitionData const& phase_transition_data,
    SpecificBodyForceData<DisplacementDim> const& specific_body_force,
    ViscosityData const& viscosity_data,
    FT2DerivativeData<DisplacementDim>& dfT_2) const
{
    auto const k_over_mu_G =
        permeability_data.Ki * permeability_data.k_rel_G / viscosity_data.mu_GR;
    auto const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    dfT_2.dp_GR_Npart = phase_transition_data.drho_GR_dp_GR *
                            fluid_enthalpy_data.h_G * darcy_velocity_data.w_GS +
                        fluid_density_data.rho_GR * fluid_enthalpy_data.h_G *
                            k_over_mu_G * phase_transition_data.drho_GR_dp_GR *
                            specific_body_force();
    dfT_2.dp_GR_gradNpart =
        fluid_density_data.rho_GR * fluid_enthalpy_data.h_G * k_over_mu_G -
        fluid_density_data.rho_LR * fluid_enthalpy_data.h_L * k_over_mu_L;

    // From p_LR = p_GR - p_cap it follows for
    // drho_LR/dp_GR = drho_LR/dp_LR * dp_LR/dp_GR
    //               = drho_LR/dp_LR * (dp_GR/dp_GR - dp_cap/dp_GR)
    //               = drho_LR/dp_LR * (1 - 0)
    double const drho_LR_dp_cap = -phase_transition_data.drho_LR_dp_LR;

    dfT_2.dp_cap_Npart =
        -drho_LR_dp_cap * fluid_enthalpy_data.h_L * darcy_velocity_data.w_LS -
        fluid_density_data.rho_LR * fluid_enthalpy_data.h_L * k_over_mu_L *
            drho_LR_dp_cap * specific_body_force();
    dfT_2.dp_cap_gradNpart =
        fluid_density_data.rho_LR * fluid_enthalpy_data.h_L * k_over_mu_L;

    dfT_2.dT = phase_transition_data.drho_GR_dT * fluid_enthalpy_data.h_G *
                   darcy_velocity_data.w_GS +
               fluid_density_data.rho_GR * phase_transition_data.dh_G_dT *
                   darcy_velocity_data.w_GS +
               phase_transition_data.drho_LR_dT * fluid_enthalpy_data.h_L *
                   darcy_velocity_data.w_LS +
               fluid_density_data.rho_LR * phase_transition_data.dh_L_dT *
                   darcy_velocity_data.w_LS;
    // TODO (naumov) + k_over_mu_G * drho_GR_dT * b + k_over_mu_L *
    // drho_LR_dT * b
}

template struct FT2Model<2>;
template struct FT2Model<3>;

template <int DisplacementDim>
void FT3Model<DisplacementDim>::eval(
    ConstituentDensityData const& constituent_density_data,
    DarcyVelocityData<DisplacementDim> const& darcy_velocity_data,
    DiffusionVelocityData<DisplacementDim> const& diffusion_velocity_data,
    FluidDensityData const& fluid_density_data,
    PhaseTransitionData const& phase_transition_data,
    SpecificBodyForceData<DisplacementDim> const& specific_body_force,
    FT3Data<DisplacementDim>& fT_3) const
{
    fT_3.N =
        (fluid_density_data.rho_GR * darcy_velocity_data.w_GS.transpose() +
         fluid_density_data.rho_LR * darcy_velocity_data.w_LS.transpose()) *
        specific_body_force();

    fT_3.gradN.noalias() =
        constituent_density_data.rho_C_GR * phase_transition_data.hCG *
            diffusion_velocity_data.d_CG +
        constituent_density_data.rho_W_GR * phase_transition_data.hWG *
            diffusion_velocity_data.d_WG;
}

template struct FT3Model<2>;
template struct FT3Model<3>;

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
