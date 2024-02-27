/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TRMVaporDiffusion.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void TRMVaporDiffusionData<DisplacementDim>::setZero()
{
    heat_capacity_vapor = 0;
    vapor_flux = GlobalDimVector<DisplacementDim>::Zero(DisplacementDim);
    storage_coefficient_by_water_vapor = 0;

    J_pT_X_dNTdN = 0;
    K_pp_X_dNTdN = 0;
    K_TT_X_dNTdN = 0;
    K_Tp_X_dNTdN = 0;
    M_Tp_X_NTN = 0;
    M_TT_X_NTN = 0;
    M_pT_X_NTN = 0;
}

template <int DisplacementDim>
void TRMVaporDiffusionModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    LiquidDensityData const& rho_L_data, SaturationData const& S_L_data,
    SaturationDataDeriv const& dS_L_data, PorosityData const& poro_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    TemperatureData<DisplacementDim> const& T_data,
    TRMVaporDiffusionData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = T_data.T;
    variables.liquid_phase_pressure = -p_cap_data.p_cap;
    // setting pG to 1 atm
    // TODO : rewrite equations s.t. p_L = pG-p_cap
    variables.gas_phase_pressure = 1.0e5;
    variables.density = rho_L_data.rho_LR;
    variables.liquid_saturation = S_L_data.S_L;

    auto const& medium = media_data.medium;

    MPL::Phase const* gas_phase =
        medium.hasPhase("Gas") ? &medium.phase("Gas") : nullptr;

    out.setZero();

    if (gas_phase && S_L_data.S_L < 1.0)
    {
        double const rho_wv =
            gas_phase->property(MaterialPropertyLib::density)
                .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

        double const drho_wv_dT =
            gas_phase->property(MaterialPropertyLib::density)
                .template dValue<double>(variables, MPL::Variable::temperature,
                                         x_t.x, x_t.t, x_t.dt);
        double const drho_wv_dp =
            gas_phase->property(MaterialPropertyLib::density)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_phase_pressure,
                                         x_t.x, x_t.t, x_t.dt);
        auto const f_Tv =
            gas_phase
                ->property(
                    MPL::PropertyType::thermal_diffusion_enhancement_factor)
                .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

        double const phi = poro_data.phi;
        variables.porosity = phi;

        double const S_g = 1.0 - S_L_data.S_L;
        double const tortuosity =
            medium.property(MPL::PropertyType::tortuosity)
                .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
        double const D_v =
            phi * S_g * tortuosity *
            gas_phase->property(MPL::PropertyType::diffusion)
                .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

        out.J_pT_X_dNTdN = f_Tv * D_v * drho_wv_dT;
        out.K_pp_X_dNTdN = D_v * drho_wv_dp;

        out.vapor_flux = -(out.J_pT_X_dNTdN * T_data.grad_T -
                           out.K_pp_X_dNTdN * p_cap_data.grad_p_cap);
        out.heat_capacity_vapor =
            gas_phase
                ->property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

        out.M_TT_X_NTN += out.heat_capacity_vapor * rho_wv * S_g * phi;

        out.storage_coefficient_by_water_vapor =
            phi * (rho_wv * dS_L_data.dS_L_dp_cap + S_g * drho_wv_dp);

        out.M_pT_X_NTN += phi * S_g * drho_wv_dT;

        //
        // Latent heat term
        //
        if (gas_phase->hasProperty(MPL::PropertyType::specific_latent_heat))
        {
            double const factor = phi * S_g / rho_L_data.rho_LR;
            // The volumetric latent heat of vaporization of liquid water
            double const L0 =
                gas_phase->property(MPL::PropertyType::specific_latent_heat)
                    .template value<double>(variables, x_t.x, x_t.t, x_t.dt) *
                rho_L_data.rho_LR;

            double const rho_wv_over_rho_L = rho_wv / rho_L_data.rho_LR;

            out.M_TT_X_NTN +=
                factor * L0 *
                (drho_wv_dT - rho_wv_over_rho_L * rho_L_data.drho_LR_dT);
            out.M_Tp_X_NTN =
                (factor * L0 *
                     (drho_wv_dp - rho_wv_over_rho_L * rho_L_data.drho_LR_dp) +
                 L0 * phi * rho_wv_over_rho_L * dS_L_data.dS_L_dp_cap);
            out.K_TT_X_dNTdN = L0 * out.J_pT_X_dNTdN / rho_L_data.rho_LR;
            out.K_Tp_X_dNTdN = L0 * out.K_pp_X_dNTdN / rho_L_data.rho_LR;
        }
    }
}

template struct TRMVaporDiffusionData<2>;
template struct TRMVaporDiffusionData<3>;
template struct TRMVaporDiffusionModel<2>;
template struct TRMVaporDiffusionModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
