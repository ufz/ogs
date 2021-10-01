/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RelPermVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermVanGenuchten::RelPermVanGenuchten(
    std::string name,
    double const residual_liquid_saturation,
    double const residual_gas_saturation,
    double const min_relative_permeability_liquid,
    double const exponent)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(1. - residual_gas_saturation),
      k_rel_min_(min_relative_permeability_liquid),
      m_(exponent)
{
    name_ = std::move(name);

    if (!(m_ > 0 && m_ < 1))
    {
        OGS_FATAL(
            "The exponent value m = {:g} of van Genuchten relative "
            "permeability model, is out of its range of (0, 1)",
            m_);
    }
}

PropertyDataType RelPermVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        S_L_res_, S_L_max_);

    double const S_eff = (S_L - S_L_res_) / (S_L_max_ - S_L_res_);
    double const v = 1. - std::pow(1. - std::pow(S_eff, 1. / m_), m_);
    double const k_rel = std::sqrt(S_eff) * v * v;
    return std::max(k_rel_min_, k_rel);
}

PropertyDataType RelPermVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (primary_variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "RelativePermeabilityVanGenuchten::dValue is implemented for "
            "derivatives with respect to liquid saturation only.");
    }

    double const S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        S_L_res_, S_L_max_);

    double const S_eff = (S_L - S_L_res_) / (S_L_max_ - S_L_res_);
    if (S_eff <= 0)  // prevent division by zero
    {
        return 0.;
    }

    if (S_eff >= 1)  // prevent taking root of zero
    {
        return 0.;
    }

    double const S_eff_to_1_over_m = std::pow(S_eff, 1. / m_);
    double const v = 1. - std::pow(1. - S_eff_to_1_over_m, m_);
    double const sqrt_S_eff = std::sqrt(S_eff);
    double const k_rel = sqrt_S_eff * v * v;

    if (k_rel < k_rel_min_)
    {
        return 0.;
    }

    return (0.5 * v * v / sqrt_S_eff +
            2. * sqrt_S_eff * v * std::pow(1. - S_eff_to_1_over_m, m_ - 1.) *
                S_eff_to_1_over_m / S_eff) /
           (S_L_max_ - S_L_res_);
}

}  // namespace MaterialPropertyLib
