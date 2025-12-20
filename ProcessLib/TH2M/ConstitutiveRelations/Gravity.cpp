// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Gravity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void GravityModel<DisplacementDim>::eval(
    FluidDensityData const& fluid_density_data,
    PorosityData const& porosity_data,
    SaturationData const& S_L_data,
    SolidDensityData const& solid_density_data,
    SpecificBodyForceData<DisplacementDim> const& specific_body_force,
    VolumetricBodyForce<DisplacementDim>& volumetric_body_force) const
{
    auto const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;
    auto const phi_L = S_L_data.S_L * porosity_data.phi;
    auto const phi_S = 1. - porosity_data.phi;

    auto const rhoGR = fluid_density_data.rho_GR;
    auto const rhoLR = fluid_density_data.rho_LR;
    auto const rho =
        phi_G * rhoGR + phi_L * rhoLR + phi_S * solid_density_data.rho_SR;
    *volumetric_body_force = rho * specific_body_force();
}

template struct GravityModel<2>;
template struct GravityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
