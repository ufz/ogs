/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EqT.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void EqTModel<DisplacementDim>::eval(
    TRMHeatStorageAndFluxData<DisplacementDim> const& heat_data,
    TRMVaporDiffusionData<DisplacementDim> const& vap_data,
    EqTData<DisplacementDim>& out) const
{
    out.M_TT_X_NTN = heat_data.M_TT_X_NTN + vap_data.M_TT_X_NTN;

    out.K_TT_NT_V_dN = heat_data.advective_heat_flux_contribution_to_K_liquid +
                       vap_data.vapor_flux * vap_data.heat_capacity_vapor;
}

template struct EqTModel<2>;
template struct EqTModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
