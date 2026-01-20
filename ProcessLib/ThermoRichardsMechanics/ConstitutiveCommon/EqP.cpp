// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EqP.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void EqPModel<DisplacementDim>::eval(
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    TemperatureData<DisplacementDim> const& T_data,
    SaturationData const& S_L_data,
    SaturationDataDeriv const& dS_L_data,
    BiotData const& biot_data,
    LiquidDensityData const& rho_L_data,
    LiquidViscosityData const& mu_L_data,
    PermeabilityData<DisplacementDim> const& perm_data,
    FluidThermalExpansionData const& f_therm_exp_data,
    TRMVaporDiffusionData<DisplacementDim> const& vap_data,
    TRMStorageData const& storage_data,
    EqPData<DisplacementDim>& out) const
{
    out.M_pu_X_BTI2N = S_L_data.S_L * rho_L_data.rho_LR * biot_data();

    out.K_pp_Laplace =
        perm_data.k_rel * rho_L_data.rho_LR * perm_data.Ki / mu_L_data();

    out.J_pp_X_BTI2NT_u_dot_N =
        -rho_L_data.rho_LR * dS_L_data.dS_L_dp_cap * biot_data();

    out.J_pp_dNT_V_N =
        perm_data.Ki / mu_L_data() *
        (rho_L_data.rho_LR * perm_data.dk_rel_dS_L * dS_L_data.dS_L_dp_cap *
         (p_cap_data.grad_p_cap + rho_L_data.rho_LR * specific_body_force_()));

    out.M_pT_X_NTN = -S_L_data.S_L * rho_L_data.rho_LR *
                         f_therm_exp_data.eff_thermal_expansion +
                     vap_data.M_pT_X_NTN;

    out.storage_p_a_p_X_NTN = storage_data.storage_p_a_p +
                              vap_data.storage_coefficient_by_water_vapor;

    out.rhs_p_dNT_V =
        -rho_L_data.rho_LR * (out.K_pp_Laplace * specific_body_force_()) +
        vap_data.J_pT_X_dNTdN * T_data.grad_T;
}

template struct EqPModel<2>;
template struct EqPModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
