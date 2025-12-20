// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DarcyLaw.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void DarcyLawModel<DisplacementDim>::eval(
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    LiquidDensityData const& rho_L_data,
    LiquidViscosityData const& mu_L_data,
    PermeabilityData<DisplacementDim> const& perm_data,
    ThermoOsmosisData<DisplacementDim> const& th_osmosis_data,
    DarcyLawData<DisplacementDim>& out) const
{
    *out = perm_data.Ki / mu_L_data() *
               (perm_data.k_rel *
                (p_cap_data.grad_p_cap + rho_L_data.rho_LR * b_)) +
           th_osmosis_data.seepage_velocity_contribution;
}

template struct DarcyLawModel<2>;
template struct DarcyLawModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
