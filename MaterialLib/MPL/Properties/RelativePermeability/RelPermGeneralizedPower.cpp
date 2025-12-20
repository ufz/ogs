// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "RelPermGeneralizedPower.h"

#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermGeneralizedPower::RelPermGeneralizedPower(
    std::string name,
    const double residual_liquid_saturation,
    const double residual_gas_saturation,
    const double min_relative_permeability,
    const double a,
    const double lambda)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      min_relative_permeability_(min_relative_permeability),
      a_(a),
      lambda_(lambda)
{
    name_ = std::move(name);
}

PropertyDataType RelPermGeneralizedPower::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = variable_array.liquid_saturation;

    if (std::isnan(S_L))
    {
        OGS_FATAL(
            "In RelPermGeneralizedPower::value, the liquid saturation is "
            "NaN.");
    }

    auto const S_L_res = residual_liquid_saturation_;
    auto const S_L_max = 1. - residual_gas_saturation_;

    auto const S_e = (S_L - S_L_res) / (S_L_max - S_L_res);

    if (S_e >= 1.0)
    {
        // fully saturated medium
        return a_;
    }
    if (S_e <= 0.0)
    {
        // dry medium
        return min_relative_permeability_;
    }

    return std::max(a_ * std::pow(S_e, lambda_), min_relative_permeability_);
}
PropertyDataType RelPermGeneralizedPower::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "RelPermGeneralizedPower::dValue is implemented for "
            "derivatives with respect to liquid saturation only.");
    }

    const double S_L = variable_array.liquid_saturation;

    auto const S_L_res = residual_liquid_saturation_;
    auto const S_L_max = 1. - residual_gas_saturation_;
    auto const S_e = (S_L - S_L_res) / (S_L_max - S_L_res);

    if ((S_e < 0.) || (S_e > 1.))
    {
        return 0.;
    }

    auto const k_rel = a_ * std::pow(S_e, lambda_);
    if (k_rel < min_relative_permeability_)
    {
        return 0.;
    }

    auto const dS_e_dS_L = 1. / (S_L_max - S_L_res);

    auto const dk_rel_dS_e = lambda_ * k_rel / S_e;
    return dk_rel_dS_e * dS_e_dS_L;
}

}  // namespace MaterialPropertyLib
