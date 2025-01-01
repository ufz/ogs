/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
    double const pressure_exponent,
    double const saturation_exponent,
    double const p_b)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(1. - residual_gas_saturation),
      m_(pressure_exponent),
      n_(saturation_exponent),
      p_b_(p_b)
{
    name_ = std::move(name);

    if (!(m_ > 0 && m_ < 1))
    {
        OGS_FATAL(
            "The pressure exponent value m = {:g} of van Genuchten saturation "
            "model, is out of its range of (0, 1)",
            m_);
    }

    if (n_ < 1)
    {
        OGS_FATAL(
            "The saturation exponent value n = {:g} of van Genuchten "
            "saturation model, is out of its range of [1, +inf)",
            n_);
    }
}

PropertyDataType SaturationVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double p_cap = variable_array.capillary_pressure;

    if (p_cap <= 0)
    {
        return S_L_max_;
    }

    double const p = p_cap / p_b_;
    double const p_to_n = std::pow(p, n_);

    double const S_eff = std::pow(p_to_n + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;
    return std::clamp(S, S_L_res_, S_L_max_);
}

PropertyDataType SaturationVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "SaturationVanGenuchten::dValue is implemented for derivatives "
            "with respect to capillary pressure only.");
    }

    const double p_cap = variable_array.capillary_pressure;

    if (p_cap <= 0)
    {
        return 0.;
    }

    double const p = p_cap / p_b_;
    double const p_to_n = std::pow(p, n_);

    double const S_eff = std::pow(p_to_n + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    if (S < S_L_res_ || S > S_L_max_)
    {
        return 0.;
    }

    double const dS_eff_dp_cap =
        -m_ * n_ * p_to_n * S_eff / (p_cap * (p_to_n + 1.));
    return dS_eff_dp_cap * (S_L_max_ - S_L_res_);
}

PropertyDataType SaturationVanGenuchten::d2Value(
    VariableArray const& variable_array, Variable const variable1,
    Variable const variable2, ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/, double const /*dt*/) const
{
    (void)variable1;
    (void)variable2;
    assert((variable1 == Variable::capillary_pressure) &&
           (variable2 == Variable::capillary_pressure) &&
           "SaturationVanGenuchten::d2Value is implemented for  derivatives "
           "with respect to capillary pressure only.");

    const double p_cap = variable_array.capillary_pressure;

    if (p_cap <= 0)
    {
        return 0.;
    }

    double const p = p_cap / p_b_;
    double const p_to_n = std::pow(p, n_);

    double const S_eff = std::pow(p_to_n + 1., -m_);
    double const S = S_eff * S_L_max_ - S_eff * S_L_res_ + S_L_res_;

    if (S < S_L_res_ || S > S_L_max_)
    {
        return 0.;
    }

    double const d2S_eff_dp_cap2 =
        m_ * n_ * n_ * p_to_n * S_eff * (p_to_n - m_) /
        ((p_cap * p_to_n + p_cap) * (p_cap * p_to_n + p_cap));
    return d2S_eff_dp_cap2 * (S_L_max_ - S_L_res_);
}
}  // namespace MaterialPropertyLib
