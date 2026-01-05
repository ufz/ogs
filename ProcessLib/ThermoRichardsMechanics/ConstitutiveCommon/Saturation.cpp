// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Saturation.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void SaturationModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    SaturationData& S_L_data, SaturationDataDeriv& dS_L_data) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.capillary_pressure = p_cap_data.p_cap;

    auto const& medium = media_data.medium;

    S_L_data.S_L = medium.property(MPL::PropertyType::saturation)
                       .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    dS_L_data.dS_L_dp_cap =
        medium.property(MPL::PropertyType::saturation)
            .template dValue<double>(variables,
                                     MPL::Variable::capillary_pressure, x_t.x,
                                     x_t.t, x_t.dt);
}

template struct SaturationModel<2>;
template struct SaturationModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
