/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Advection.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void AdvectionModel<DisplacementDim>::eval(
    ConstituentDensityData const& constituent_density_data,
    PermeabilityData<DisplacementDim> const& permeability_data,
    PureLiquidDensityData const& rho_W_LR,
    ViscosityData const& viscosity_data,
    AdvectionData<DisplacementDim>& advection_data) const
{
    GlobalDimMatrix<DisplacementDim> const k_over_mu_G =
        permeability_data.Ki * permeability_data.k_rel_G / viscosity_data.mu_GR;
    GlobalDimMatrix<DisplacementDim> const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    advection_data.advection_C_G =
        constituent_density_data.rho_C_GR * k_over_mu_G;
    advection_data.advection_C_L =
        constituent_density_data.rho_C_LR * k_over_mu_L;
    advection_data.advection_W_G =
        constituent_density_data.rho_W_GR * k_over_mu_G;
    advection_data.advection_W_L = rho_W_LR() * k_over_mu_L;
}

template <int DisplacementDim>
void AdvectionModel<DisplacementDim>::dEval(
    ConstituentDensityData const& constituent_density_data,
    PermeabilityData<DisplacementDim> const& permeability_data,
    ViscosityData const& viscosity_data,
    SaturationDataDeriv const& dS_L_dp_cap,
    PhaseTransitionData const& phase_transition_data,
    AdvectionDerivativeData<DisplacementDim>& advection_d_data) const
{
    GlobalDimMatrix<DisplacementDim> const k_over_mu_G =
        permeability_data.Ki * permeability_data.k_rel_G / viscosity_data.mu_GR;
    GlobalDimMatrix<DisplacementDim> const k_over_mu_L =
        permeability_data.Ki * permeability_data.k_rel_L / viscosity_data.mu_LR;

    GlobalDimMatrix<DisplacementDim> const dk_over_mu_G_dp_cap =
        permeability_data.Ki * permeability_data.dk_rel_G_dS_L * dS_L_dp_cap() /
        viscosity_data.mu_GR;
    GlobalDimMatrix<DisplacementDim> const dk_over_mu_L_dp_cap =
        permeability_data.Ki * permeability_data.dk_rel_L_dS_L * dS_L_dp_cap() /
        viscosity_data.mu_LR;

    advection_d_data.dadvection_C_dp_GR =
        phase_transition_data.drho_C_GR_dp_GR * k_over_mu_G
        // + rhoCGR * (dk_over_mu_G_dp_GR = 0)
        // + rhoCLR * (dk_over_mu_L_dp_GR = 0)
        + phase_transition_data.drho_C_LR_dp_GR * k_over_mu_L;

    advection_d_data.dadvection_C_dp_cap =
        //(drho_C_GR_dp_cap = 0) * k_over_mu_G
        constituent_density_data.rho_C_GR * dk_over_mu_G_dp_cap +
        (-phase_transition_data.drho_C_LR_dp_LR) * k_over_mu_L +
        constituent_density_data.rho_C_LR * dk_over_mu_L_dp_cap;
}

template struct AdvectionModel<2>;
template struct AdvectionModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
