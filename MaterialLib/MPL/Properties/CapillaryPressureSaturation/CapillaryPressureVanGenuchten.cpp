/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 20, 2020, 9:59 AM
 */

#include "CapillaryPressureVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
CapillaryPressureVanGenuchten::CapillaryPressureVanGenuchten(
    std::string name,
    double const residual_liquid_saturation,
    double const residual_gas_saturation,
    double const exponent,
    double const p_b,
    double const maximum_capillary_pressure)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(1. - residual_gas_saturation),
      m_(exponent),
      p_b_(p_b),
      p_cap_max_(maximum_capillary_pressure)
{
    name_ = std::move(name);

    if (S_L_res_ < 0 || S_L_res_ > 1)
    {
        OGS_FATAL(
            "Van Genuchten capillary pressure model: "
            "The residual liquid saturation value S_L_res = {:g} is out of "
            "admissible range [0, 1].",
            S_L_res_);
    }
    if (S_L_max_ < 0 || S_L_max_ > 1)
    {
        OGS_FATAL(
            "Van Genuchten capillary pressure model: "
            "The maximum liquid saturation value S_L_max = {:g} is out of "
            "admissible range [0, 1].",
            S_L_max_);
    }
    if (S_L_res_ >= S_L_max_)
    {
        OGS_FATAL(
            "Van Genuchten capillary pressure model: "
            "The maximum liquid saturation S_L_max = {:g} must not be less or "
            "equal to the residual liquid saturion S_L_res = { : g}.",
            S_L_max_, S_L_res_);
    }
    if (!(m_ > 0 && m_ < 1))
    {
        OGS_FATAL(
            "Van Genuchten capillary pressure model: "
            "The exponent value m = {:g} is out of of admissible range (0, 1).",
            m_);
    }
    if (p_b_ <= 0)
    {
        OGS_FATAL(
            "Van Genuchten capillary pressure model: "
            "The pressure scaling value p_b = {:g} must be positive.",
            p_b_);
    }
    if (p_cap_max_ < 0)
    {
        OGS_FATAL(
            "Van Genuchten capillary pressure model: "
            "The maximum capillary pressure value p_cap_max = {:g} must be "
            "non-negative.",
            p_cap_max_);
    }
}

PropertyDataType CapillaryPressureVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (S_L <= S_L_res_)
    {
        return p_cap_max_;
    }

    if (S_L >= S_L_max_)
    {
        return 0.;
    }

    double const S_eff = (S_L - S_L_res_) / (S_L_max_ - S_L_res_);
    assert(0 <= S_eff && S_eff <= 1);

    double const p_cap =
        p_b_ * std::pow(std::pow(S_eff, -1.0 / m_) - 1.0, 1.0 - m_);
    assert(p_cap > 0);
    return std::min(p_cap, p_cap_max_);
}

PropertyDataType CapillaryPressureVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "CapillaryPressureVanGenuchten::dValue is implemented for "
           "derivatives with respect to liquid saturation only.");

    double const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (S_L <= S_L_res_)
    {
        return 0.;
    }

    if (S_L >= S_L_max_)
    {
        return 0.;
    }

    double const S_eff = (S_L - S_L_res_) / (S_L_max_ - S_L_res_);

    assert(0 < S_eff && S_eff < 1.0);

    double const val1 = std::pow(S_eff, -1.0 / m_);
    double const p_cap = p_b_ * std::pow(val1 - 1.0, 1.0 - m_);
    if (p_cap >= p_cap_max_)
    {
        return 0.;
    }

    double const val2 = std::pow(val1 - 1.0, -m_);
    return p_b_ * (m_ - 1.0) * val1 * val2 / (m_ * (S_L - S_L_res_));
}
}  // namespace MaterialPropertyLib
