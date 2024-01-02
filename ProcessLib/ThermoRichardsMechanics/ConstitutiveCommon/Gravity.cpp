/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Gravity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void GravityModel<DisplacementDim>::eval(
    PorosityData const& poro_data, SolidDensityData const& rho_S_data,
    LiquidDensityData const& rho_L_data, SaturationData const& S_L_data,
    SaturationDataDeriv const& dS_L_data,
    GravityData<DisplacementDim>& out) const
{
    auto const rho_SR = rho_S_data.rho_SR;
    auto const phi = poro_data.phi;
    auto const S_L = S_L_data.S_L;
    auto const rho_LR = rho_L_data.rho_LR;
    auto const b = specific_body_force_;

    double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
    out.volumetric_body_force = rho * b;

    // There is no minus in the Jacobian block, because volumetric_body_force is
    // subtracted from the residual vector.
    out.J_up_HT_V_N =
        poro_data.phi * rho_L_data.rho_LR * dS_L_data.dS_L_dp_cap * b;
}

template struct GravityModel<2>;
template struct GravityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
