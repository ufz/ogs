// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "ConstitutiveData.h"
#include "ConstitutiveModels.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/MaterialState.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
struct ConstitutiveSetting
{
    /// Evaluate the constitutive setting.
    void eval(
        ConstitutiveModels<DisplacementDim>& models, double const t,
        double const dt, ParameterLib::SpatialPosition const& x_position,
        MaterialPropertyLib::Medium const& medium, TemperatureData const T_data,
        CapillaryPressureData<DisplacementDim> const& p_cap_data,
        KelvinVector<DisplacementDim> const& eps_arg,
        StatefulData<DisplacementDim>& state,
        StatefulDataPrev<DisplacementDim> const& prev_state,
        ProcessLib::ThermoRichardsMechanics::MaterialStateData<DisplacementDim>&
            mat_state,
        ConstitutiveTempData<DisplacementDim>& tmp,
        OutputData<DisplacementDim>& out,
        ConstitutiveData<DisplacementDim>& cd) const;
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ProcessLib::RichardsMechanics
