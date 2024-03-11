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
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
struct ConstitutiveSetting
{
    void init(ConstitutiveModels<DisplacementDim>&, double const /*t*/,
              double const /*dt*/, ParameterLib::SpatialPosition const&,
              MediaData const&, TemperatureData<DisplacementDim> const&,
              StatefulData<DisplacementDim>&,
              StatefulDataPrev<DisplacementDim>&) const
    {
    }

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
        return std::get<TotalStressData<DisplacementDim>>(state).sigma_total;
    }
    static KelvinVector<DisplacementDim>& statefulStress(
        StatefulData<DisplacementDim>& state)
    {
        return std::get<TotalStressData<DisplacementDim>>(state).sigma_total;
    }

    /// In case that the input initial data for
    /// state.s_mech_data.sigma_total are effective stress values,
    /// state.s_mech_data.sigma_total is reset to total stress.
    static void convertInitialStressType(
        StatefulData<DisplacementDim>& state,
        StatefulDataPrev<DisplacementDim>& prev_state,
        KelvinVector<DisplacementDim> const& pore_pressure_part)
    {
        auto& sigma_total =
            std::get<TotalStressData<DisplacementDim>>(state).sigma_total;
        sigma_total.noalias() -= pore_pressure_part;

        (std::get<PrevState<TotalStressData<DisplacementDim>>>(prev_state)
             ->sigma_total)
            .noalias() = sigma_total;
    }
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
