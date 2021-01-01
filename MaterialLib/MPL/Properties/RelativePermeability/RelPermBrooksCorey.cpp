/**
 * \file
 * \author Norbert Grunwald
 * \date   02.07.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RelPermBrooksCorey.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermBrooksCorey::RelPermBrooksCorey(
    std::string name,
    const double residual_liquid_saturation,
    const double residual_gas_saturation,
    const double min_relative_permeability_liquid,
    const double min_relative_permeability_gas,
    const double exponent)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      min_relative_permeability_liquid_(min_relative_permeability_liquid),
      min_relative_permeability_gas_(min_relative_permeability_gas),
      exponent_(exponent)
{
    name_ = std::move(name);
};

PropertyDataType RelPermBrooksCorey::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    /// here, an extra computation of saturation is forced, guaranteeing a
    /// correct value. In order to speed up the computing time, saturation could
    /// be insertred into the primary variable array after it is computed in the
    /// FEM assembly.
    auto const s_L = std::visit(
        [&variable_array, &pos, t, dt](auto&& scale) -> double {
            return scale->property(PropertyType::saturation)
                .template value<double>(variable_array, pos, t, dt);
        },
        scale_);

    auto const s_L_res = residual_liquid_saturation_;
    auto const s_L_max = 1. - residual_gas_saturation_;
    auto const k_rel_min_LR = min_relative_permeability_liquid_;
    auto const k_rel_min_GR = min_relative_permeability_gas_;

    auto const lambda = exponent_;

    auto const s_eff = (s_L - s_L_res) / (s_L_max - s_L_res);

    if (s_eff >= 1.0)
    {
        // fully saturated medium
        return Eigen::Vector2d{1.0, k_rel_min_GR};
    }
    if (s_eff <= 0.0)
    {
        // dry medium
        return Eigen::Vector2d{k_rel_min_LR, 1.0};
    }

    auto const k_rel_LR = std::pow(s_eff, (2. + 3. * lambda) / lambda);
    auto const k_rel_GR = (1. - s_eff) * (1. - s_eff) *
                          (1. - std::pow(s_eff, (2. + lambda) / lambda));

    return Eigen::Vector2d{std::max(k_rel_LR, k_rel_min_LR),
                           std::max(k_rel_GR, k_rel_min_GR)};
}
PropertyDataType RelPermBrooksCorey::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "RelPermBrooksCorey::dValue is implemented for "
           " derivatives with respect to liquid saturation only.");
    /// here, an extra computation of saturation is forced, guaranteeing a
    /// correct value. In order to speed up the computing time, saturation could
    /// be insertred into the primary variable array after it is computed in the
    /// FEM assembly.
    auto const s_L = std::visit(
        [&variable_array, &pos, t, dt](auto&& scale) -> double {
            return scale->property(PropertyType::saturation)
                .template value<double>(variable_array, pos, t, dt);
        },
        scale_);

    auto const s_L_res = residual_liquid_saturation_;
    auto const s_L_max = 1. - residual_gas_saturation_;
    auto const lambda = exponent_;

    auto const s_eff = (s_L - s_L_res) / (s_L_max - s_L_res);
    if ((s_eff < 0.) || (s_eff > 1.))
        return Eigen::Vector2d{0., 0.};

    auto const d_se_d_sL = 1. / (s_L_max - s_L_res);
    auto const dk_rel_LRdse =
        (3 * lambda + 2.) / lambda * std::pow(s_eff, 2. / lambda + 2.);

    auto const dk_rel_LRdsL = dk_rel_LRdse * d_se_d_sL;

    auto const twoL_L = (2. + lambda) / lambda;
    auto const dk_rel_GRdse =
        -2. * (1 - s_eff) * (1. - std::pow(s_eff, twoL_L)) -
        twoL_L * std::pow(s_eff, twoL_L - 1.) * (1. - s_eff) * (1. - s_eff);

    auto const dk_rel_GRdsL = dk_rel_GRdse * d_se_d_sL;
    return Eigen::Vector2d{dk_rel_LRdsL, dk_rel_GRdsL};
}

}  // namespace MaterialPropertyLib
