// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ConstitutiveSetting.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::eval(
    ConstitutiveModels<DisplacementDim>& /*models*/, double const /*t*/,
    double const /*dt*/, ParameterLib::SpatialPosition const& /*x_position*/,
    MaterialPropertyLib::Medium const& /*medium*/,
    TemperatureData const /*T_data*/,
    CapillaryPressureData<DisplacementDim> const& /*p_cap_data*/,
    KelvinVector<DisplacementDim> const& /*eps_arg*/,
    StatefulData<DisplacementDim>& /*state*/,
    StatefulDataPrev<DisplacementDim> const& /*prev_state*/,
    ProcessLib::ThermoRichardsMechanics::MaterialStateData<DisplacementDim>&
    /*mat_state*/,
    ConstitutiveTempData<DisplacementDim>& /*tmp*/,
    OutputData<DisplacementDim>& /*out*/,
    ConstitutiveData<DisplacementDim>& /*cd*/) const
{
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ProcessLib::RichardsMechanics
