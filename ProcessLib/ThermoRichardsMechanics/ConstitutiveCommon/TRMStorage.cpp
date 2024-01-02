/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TRMStorage.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void TRMStorageModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, BiotData const& biot_data,
    PorosityData const& poro_data, LiquidDensityData const& rho_L_data,
    SaturationData const& S_L_data, SaturationDataDeriv const& dS_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    SolidCompressibilityData const& solid_compressibility_data,
    TRMStorageData& out) const
{
    double const p_cap = p_cap_data.p_cap;
    double const p_cap_prev = p_cap_data.p_cap_prev;

    double const phi = poro_data.phi;
    double const alphaB_minus_phi = biot_data() - phi;

    double const beta_LR = rho_L_data.drho_LR_dp / rho_L_data.rho_LR;

    double const a0 = alphaB_minus_phi * solid_compressibility_data.beta_SR;
    double const specific_storage_a_p =
        S_L_data.S_L * (phi * beta_LR + S_L_data.S_L * a0);
    double const specific_storage_a_S = phi - p_cap * S_L_data.S_L * a0;

    // Note: d beta_LR/d p is omitted because it is a small value.
    double const dspecific_storage_a_p_dp_cap =
        dS_L_data.dS_L_dp_cap * (phi * beta_LR + 2 * S_L_data.S_L * a0);
    double const dspecific_storage_a_S_dp_cap =
        -a0 * (S_L_data.S_L + p_cap * dS_L_data.dS_L_dp_cap);

    // secant derivative from time discretization for storage
    // use tangent, if secant is not available
    double const DeltaS_L_Deltap_cap =
        (p_cap == p_cap_prev)
            ? dS_L_data.dS_L_dp_cap
            : (S_L_data.S_L - S_L_prev_data->S_L) / (p_cap - p_cap_prev);

    out.storage_p_a_p = rho_L_data.rho_LR * specific_storage_a_p;
    out.storage_p_a_S_X_NTN =
        -rho_L_data.rho_LR * specific_storage_a_S * DeltaS_L_Deltap_cap;
    out.J_pp_X_NTN = (p_cap - p_cap_prev) / x_t.dt * rho_L_data.rho_LR *
                     dspecific_storage_a_p_dp_cap;
    out.storage_p_a_S_Jpp_X_NTN =
        -rho_L_data.rho_LR *
        ((S_L_data.S_L - S_L_prev_data->S_L) * dspecific_storage_a_S_dp_cap +
         specific_storage_a_S * dS_L_data.dS_L_dp_cap) /
        x_t.dt;
}

template struct TRMStorageModel<2>;
template struct TRMStorageModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
