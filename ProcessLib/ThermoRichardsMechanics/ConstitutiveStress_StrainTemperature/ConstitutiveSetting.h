/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
struct ConstitutiveSetting
{
    void init(ConstitutiveModels<DisplacementDim>& models, double const t,
              double const dt, ParameterLib::SpatialPosition const& x_position,
              MediaData const& media_data,
              TemperatureData<DisplacementDim> const& T_data,
              StatefulData<DisplacementDim>& state,
              StatefulDataPrev<DisplacementDim>& prev_state) const;

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

    static KelvinVector<DisplacementDim> const& statefulStress(
        StatefulData<DisplacementDim> const& state)
    {
        return std::get<EffectiveStressData<DisplacementDim>>(state).sigma_eff;
    }
    static KelvinVector<DisplacementDim>& statefulStress(
        StatefulData<DisplacementDim>& state)
    {
        return std::get<EffectiveStressData<DisplacementDim>>(state).sigma_eff;
    }

    /// In case that the input initial data for
    /// state.s_mech_data.sigma_eff are total stress values,
    /// state.s_mech_data.sigma_eff is reset to effective stress.
    static void convertInitialStressType(
        StatefulData<DisplacementDim>& state,
        StatefulDataPrev<DisplacementDim>& prev_state,
        KelvinVector<DisplacementDim> const& pore_pressure_part)
    {
        auto& sigma_eff =
            std::get<EffectiveStressData<DisplacementDim>>(state).sigma_eff;
        sigma_eff.noalias() += pore_pressure_part;

        (std::get<PrevState<EffectiveStressData<DisplacementDim>>>(prev_state)
             ->sigma_eff)
            .noalias() = sigma_eff;
    }
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
