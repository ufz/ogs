// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DarcyVelocity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void DarcyVelocityModel<DisplacementDim>::eval(
    CapillaryPressureGradientData<DisplacementDim> const& grad_p_cap,
    FluidDensityData const& fluid_density_data,
    GasPressureGradientData<DisplacementDim> const& grad_p_GR,
    PermeabilityData<DisplacementDim> const& permeability_data,
    SpecificBodyForce<DisplacementDim> const& specific_body_force,
    ViscosityData const& viscosity_data,
    DarcyVelocityData<DisplacementDim>& darcy_velocity_data) const
{
    auto const k_over_mu_G =
        permeability_data.Ki * permeability_data.k_rel_G / viscosity_data.mu_GR;
    auto const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    darcy_velocity_data.w_GS =
        k_over_mu_G * fluid_density_data.rho_GR * specific_body_force() -
        k_over_mu_G * grad_p_GR();
    darcy_velocity_data.w_LS =
        k_over_mu_L * grad_p_cap() +
        k_over_mu_L * fluid_density_data.rho_LR * specific_body_force() -
        k_over_mu_L * grad_p_GR();
}

template struct DarcyVelocityModel<2>;
template struct DarcyVelocityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
