// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DiffusionVelocity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void DiffusionVelocityModel<DisplacementDim>::eval(
    CapillaryPressureGradientData<DisplacementDim> const& grad_p_cap,
    GasPressureGradientData<DisplacementDim> const& grad_p_GR,
    MassMoleFractionsData const& mass_mole_fractions_data,
    PhaseTransitionData const& phase_transition_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    TemperatureGradientData<DisplacementDim> const& grad_T,
    DiffusionVelocityData<DisplacementDim>& diffusion_velocity_data) const
{
    auto const gradxmWG = phase_transition_data.dxmWG_dpGR * grad_p_GR() +
                          phase_transition_data.dxmWG_dpCap * grad_p_cap() +
                          phase_transition_data.dxmWG_dT * grad_T();
    auto const gradxmCG = -gradxmWG;

    auto const gradxmWL = phase_transition_data.dxmWL_dpGR * grad_p_GR() +
                          phase_transition_data.dxmWL_dpCap * grad_p_cap() +
                          phase_transition_data.dxmWL_dT * grad_T();
    auto const gradxmCL = -gradxmWL;

    double const phi_L = S_L_data.S_L * porosity_data.phi;
    double const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;

    double const phi_G_D_vapour =
        phi_G * phase_transition_data.diffusion_coefficient_vapour;
    double const phi_L_D_solute =
        phi_L * phase_transition_data.diffusion_coefficient_solute;

    if (mass_mole_fractions_data.xmCG == 0)
    {
        diffusion_velocity_data.d_CG = GlobalDimVector<DisplacementDim>::Zero();
    }
    else
    {
        diffusion_velocity_data.d_CG =
            -phi_G_D_vapour / mass_mole_fractions_data.xmCG * gradxmCG;
    }

    if (mass_mole_fractions_data.xmCG == 1)
    {
        diffusion_velocity_data.d_WG = GlobalDimVector<DisplacementDim>::Zero();
    }
    else
    {
        diffusion_velocity_data.d_WG =
            -phi_G_D_vapour / (1 - mass_mole_fractions_data.xmCG) * gradxmWG;
    }

    if (mass_mole_fractions_data.xmWL == 1)

    {
        diffusion_velocity_data.d_CL = GlobalDimVector<DisplacementDim>::Zero();
    }
    else
    {
        diffusion_velocity_data.d_CL =
            -phi_L_D_solute / (1. - mass_mole_fractions_data.xmWL) * gradxmCL;
    }

    if (mass_mole_fractions_data.xmWL == 0)
    {
        diffusion_velocity_data.d_WL = GlobalDimVector<DisplacementDim>::Zero();
    }
    else
    {
        diffusion_velocity_data.d_WL =
            -phi_L_D_solute / mass_mole_fractions_data.xmWL * gradxmWL;
    }
}

template struct DiffusionVelocityModel<2>;
template struct DiffusionVelocityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
