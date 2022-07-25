/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EqU.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void EqUModel<DisplacementDim>::eval(
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    SaturationDataDeriv const& dS_L_data,
    BiotData const& biot_data,
    BishopsData const& bishops_data,
    LiquidDensityData const& rho_L_data,
    PorosityData const& poro_data,
    EqUData<DisplacementDim>& out) const
{
    out.J_up_X_BTI2N =
        -biot_data.alpha *
        (bishops_data.chi_S_L +
         bishops_data.dchi_dS_L * p_cap_data.p_cap * dS_L_data.dS_L_dp_cap);

    out.J_up_HT_V_N =
        poro_data.phi * rho_L_data.rho_LR * dS_L_data.dS_L_dp_cap * b_;
}

template struct EqUModel<2>;
template struct EqUModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
