/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaturationVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
SaturationVanGenuchten::SaturationVanGenuchten(
    std::string name,
    double const residual_liquid_saturation,
    double const residual_gas_saturation,
    double const exponent,
    double const p_b)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(1. - residual_gas_saturation),
      m_(exponent),
      p_b_(p_b)
{
    name_ = std::move(name);

    if (!(m_ > 0 && m_ < 1))
    {
        OGS_FATAL(
            "The exponent value m = {:g} of van Genuchten saturation model, is "
            "out of its range of (0, 1)",
            m_);
    }
}

PropertyDataType SaturationVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if (p_cap <= 0)
    {
        return S_L_max_;
    }

    double const p = p_cap / p_b_;
    double const n = 1. / (1. - m_);
    double const p_to_n = std::pow(p, n);

    double const S_eff = std::pow(p_to_n + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;
    return std::clamp(S, S_L_res_, S_L_max_);
}

PropertyDataType SaturationVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (primary_variable != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "SaturationVanGenuchten::dValue is implemented for derivatives "
            "with respect to capillary pressure only.");
    }

    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if (p_cap <= 0)
    {
        return 0.;
    }

    double const p = p_cap / p_b_;
    double const n = 1. / (1. - m_);
    double const p_to_n = std::pow(p, n);

    double const S_eff = std::pow(p_to_n + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    if (S < S_L_res_ || S > S_L_max_)
    {
        return 0.;
    }

    double const dS_eff_dp_cap = -m_ * std::pow(p, n - 1) *
                                 std::pow(1 + p_to_n, -1. - m_) /
                                 (p_b_ * (1. - m_));
    return dS_eff_dp_cap * (S_L_max_ - S_L_res_);
}

PropertyDataType SaturationVanGenuchten::d2Value(
    VariableArray const& variable_array, Variable const primary_variable1,
    Variable const primary_variable2,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable1;
    (void)primary_variable2;
    assert((primary_variable1 == Variable::capillary_pressure) &&
           (primary_variable2 == Variable::capillary_pressure) &&
           "SaturationVanGenuchten::d2Value is implemented for  derivatives "
           "with respect to capillary pressure only.");

    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if (p_cap <= 0)
    {
        return 0.;
    }

    double const p = p_cap / p_b_;
    double const n = 1. / (1. - m_);
    double const p_to_n = std::pow(p, n);

    double const S_eff = std::pow(p_to_n + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    if (S < S_L_res_ || S > S_L_max_)
    {
        return 0.;
    }

    double const d2S_eff_dp_cap2 =
        m_ * p_to_n * std::pow(p_to_n + 1., -m_ - 2.) * (p_to_n - m_) /
        (p_cap * p_cap * (m_ - 1.) * (m_ - 1.));
    return d2S_eff_dp_cap2 * (S_L_max_ - S_L_res_);
}
}  // namespace MaterialPropertyLib
