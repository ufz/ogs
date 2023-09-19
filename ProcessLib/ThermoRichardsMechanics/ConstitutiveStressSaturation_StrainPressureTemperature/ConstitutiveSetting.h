/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ConstitutiveData.h"
#include "ConstitutiveModels.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
struct ConstitutiveSetting
{
    /// Evaluate the constitutive setting.
    void eval(ConstitutiveModels<DisplacementDim>& models, double const t,
              double const dt, ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium const& medium,
              TemperatureData<DisplacementDim> const& T_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              KelvinVector<DisplacementDim> const& eps_arg,
              StatefulData<DisplacementDim>& state,
              StatefulDataPrev<DisplacementDim> const& prev_state,
              MaterialStateData<DisplacementDim>& mat_state,
              ConstitutiveTempData<DisplacementDim>& tmp,
              OutputData<DisplacementDim>& out,
              ConstitutiveData<DisplacementDim>& cd) const;

    static KelvinVector<DisplacementDim> const& totalStress(
        ConstitutiveData<DisplacementDim> const& /*cd*/,
        StatefulData<DisplacementDim> const& state)
    {
        return state.total_stress_data.sigma_total;
    }
    static KelvinVector<DisplacementDim>& totalStress(
        ConstitutiveData<DisplacementDim>& /*cd*/,
        StatefulData<DisplacementDim>& state)
    {
        return state.total_stress_data.sigma_total;
    }

    static KelvinVector<DisplacementDim> const& statefulStress(
        StatefulData<DisplacementDim> const& state)
    {
        return state.total_stress_data.sigma_total;
    }
    static KelvinVector<DisplacementDim>& statefulStress(
        StatefulData<DisplacementDim>& state)
    {
        return state.total_stress_data.sigma_total;
    }
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
