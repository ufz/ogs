// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "RelPermBrooksCorey.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermBrooksCorey::RelPermBrooksCorey(std::string name,
                                       const double residual_liquid_saturation,
                                       const double residual_gas_saturation,
                                       const double min_relative_permeability,
                                       const double exponent)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      min_relative_permeability_(min_relative_permeability),
      exponent_(exponent)
{
    name_ = std::move(name);

    if (exponent_ <= 0.)
    {
        OGS_FATAL(
            "RelPermBrooksCorey: exponent 'lambda' must be positive, but {} "
            "was given.",
            exponent_);
    }
    if (residual_liquid_saturation_ < 0.)
    {
        OGS_FATAL(
            "RelPermBrooksCorey: residual_liquid_saturation must be "
            "non-negative, but {} was given.",
            residual_liquid_saturation_);
    }
    if (residual_gas_saturation_ < 0.)
    {
        OGS_FATAL(
            "RelPermBrooksCorey: residual_gas_saturation must be non-negative, "
            "but {} was given.",
            residual_gas_saturation_);
    }
    if (residual_liquid_saturation_ + residual_gas_saturation_ >= 1.)
    {
        OGS_FATAL(
            "RelPermBrooksCorey: residual_liquid_saturation ({}) + "
            "residual_gas_saturation ({}) must be less than 1 so that the "
            "effective saturation range is positive.",
            residual_liquid_saturation_, residual_gas_saturation_);
    }
    if (min_relative_permeability_ < 0. || min_relative_permeability_ > 1.)
    {
        OGS_FATAL(
            "RelPermBrooksCorey: min_relative_permeability must be in [0, 1], "
            "but {} was given.",
            min_relative_permeability_);
    }
};

PropertyDataType RelPermBrooksCorey::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    /// here, an extra computation of saturation is forced, guaranteeing a
    /// correct value. In order to speed up the computing time, saturation could
    /// be inserted into the primary variable array after it is computed in the
    /// FEM assembly.
    auto const s_L = std::visit(
        [&variable_array, &pos, t, dt](auto&& scale) -> double
        {
            return scale->property(PropertyType::saturation)
                .template value<double>(variable_array, pos, t, dt);
        },
        scale_);

    auto const s_L_res = residual_liquid_saturation_;
    auto const s_L_max = 1. - residual_gas_saturation_;

    auto const lambda = exponent_;

    auto const s_eff = (s_L - s_L_res) / (s_L_max - s_L_res);

    if (s_eff >= 1.0)
    {
        // fully saturated medium
        return 1.0;
    }
    if (s_eff <= 0.0)
    {
        // dry medium
        return min_relative_permeability_;
    }

    auto const k_rel_LR = std::pow(s_eff, (2. + 3. * lambda) / lambda);

    return std::max(k_rel_LR, min_relative_permeability_);
}
PropertyDataType RelPermBrooksCorey::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "RelPermBrooksCorey::dValue is implemented for derivatives with "
            "respect to liquid saturation only.");
    }

    /// here, an extra computation of saturation is forced, guaranteeing a
    /// correct value. In order to speed up the computing time, saturation could
    /// be inserted into the primary variable array after it is computed in the
    /// FEM assembly.
    auto const s_L = std::visit(
        [&variable_array, &pos, t, dt](auto&& scale) -> double
        {
            return scale->property(PropertyType::saturation)
                .template value<double>(variable_array, pos, t, dt);
        },
        scale_);

    auto const s_L_res = residual_liquid_saturation_;
    auto const s_L_max = 1. - residual_gas_saturation_;
    auto const lambda = exponent_;

    auto const s_eff = (s_L - s_L_res) / (s_L_max - s_L_res);
    // Consistency with value(): at s_eff == 1 the value is clamped to the
    // constant 1.0 (and saturation models clamping S_L at s_L_max put whole
    // regions exactly on this point), so the derivative there is zero.
    if ((s_eff < 0.) || (s_eff >= 1.))
    {
        return 0.;
    }

    // Consistency with value(): where the min_relative_permeability clamp is
    // active (dry range), the relative permeability is constant and its
    // derivative is zero.
    auto const k_rel_LR = std::pow(s_eff, (2. + 3. * lambda) / lambda);
    if (k_rel_LR <= min_relative_permeability_)
    {
        return 0.;
    }

    auto const d_se_d_sL = 1. / (s_L_max - s_L_res);
    auto const dk_rel_LRdse =
        (3 * lambda + 2.) / lambda * std::pow(s_eff, 2. / lambda + 2.);

    return dk_rel_LRdse * d_se_d_sL;
}

}  // namespace MaterialPropertyLib
