/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TRMHeatStorageAndFlux.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void TRMHeatStorageAndFluxModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    LiquidDensityData const& rho_L_data, SolidDensityData const& rho_S_data,
    SaturationData const& S_L_data, SaturationDataDeriv const& dS_L_data,
    PorosityData const& poro_data, LiquidViscosityData const& mu_L_data,
    PermeabilityData<DisplacementDim> const& perm,
    TemperatureData<DisplacementDim> const& T_data,
    DarcyLawData<DisplacementDim> const& darcy_data,
    TRMHeatStorageAndFluxData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = T_data.T;
    variables.porosity = poro_data.phi;
    variables.liquid_saturation = S_L_data.S_L;

    auto const& liquid_phase = media_data.liquid;
    auto const& solid_phase = media_data.solid;

    auto const specific_heat_capacity_fluid =
        liquid_phase.property(MaterialPropertyLib::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    auto const specific_heat_capacity_solid =
        solid_phase
            .property(MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    double const phi = poro_data.phi;

    // TODO real vs. "non-real" values
    double const volumetric_heat_capacity_liquid =
        rho_L_data.rho_LR * specific_heat_capacity_fluid;

    // NOTE: Gas phase is not included.
    double const volumetric_heat_capacity_liquid_and_solid =
        rho_S_data.rho_SR * specific_heat_capacity_solid * (1 - phi) +
        S_L_data.S_L * phi * volumetric_heat_capacity_liquid;

    out.M_TT_X_NTN = volumetric_heat_capacity_liquid_and_solid;

    out.K_TT_Laplace = MaterialPropertyLib::formEigenTensor<DisplacementDim>(
        media_data.medium
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variables, x_t.x, x_t.t, x_t.dt));

    // Unit is J / m^2 / s / K. It's not a heat flux, but related.
    out.advective_heat_flux_contribution_to_K_liquid =
        volumetric_heat_capacity_liquid * darcy_data.v_darcy;

    //
    // temperature equation, pressure part
    //
    out.K_Tp_NT_V_dN = -volumetric_heat_capacity_liquid * perm.k_rel /
                       mu_L_data.viscosity *
                       (perm.Ki.transpose() * T_data.grad_T);
    out.K_Tp_X_NTN = -volumetric_heat_capacity_liquid *
                     darcy_data.v_darcy.dot(T_data.grad_T) / perm.k_rel *
                     perm.dk_rel_dS_L * dS_L_data.dS_L_dp_cap;
}

template struct TRMHeatStorageAndFluxModel<2>;
template struct TRMHeatStorageAndFluxModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
