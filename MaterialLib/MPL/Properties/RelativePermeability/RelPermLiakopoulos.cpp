/**
 * \file
 * \author Norbert Grunwald
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RelPermLiakopoulos.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
/**
 */
PropertyDataType RelPermLiakopoulos::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    /// here, an extra computation of saturation is forced, guaranteeing a
    /// correct value. In order to speed up the computing time, saturation could
    /// be insertred into the primary variable array after it is computed in the
    /// FEM assembly.
    auto const s_L = _medium->property(PropertyType::saturation)
                         .template value<double>(variable_array, pos, t, dt);
    auto const s_L_res = _residual_liquid_saturation;
    auto const k_rel_min_GR = _min_relative_permeability_gas;

    if (s_L <= s_L_res)
    {
        return Eigen::Vector2d{0., 1.};
    }

    if (s_L >= 1.)
    {
        return Eigen::Vector2d{1., k_rel_min_GR};
    }

    auto const a = _parameter_a;
    auto const b = _parameter_b;
    auto const lambda = _exponent;

    auto const s_eff = (s_L - s_L_res) / (1. - s_L_res);

    auto const k_rel_LR = 1. - a * std::pow(1. - s_L, b);
    auto const k_rel_GR = (1. - s_eff) * (1. - s_eff) *
                          (1. - std::pow(s_eff, (2. + lambda) / lambda));

    return Eigen::Vector2d{std::max(k_rel_LR, 0.),
                           std::max(k_rel_GR, k_rel_min_GR)};
}

PropertyDataType RelPermLiakopoulos::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "RelPermLiakopoulos::dValue is implemented for "
           " derivatives with respect to liquid saturation only.");
    auto const s_L = _medium->property(PropertyType::saturation)
                         .template value<double>(variable_array, pos, t, dt);
    auto const s_L_res = _residual_liquid_saturation;
    auto const s_L_max = _maximal_liquid_saturation;

    const double s_L_within_range = std::min(std::max(s_L_res, s_L), s_L_max);

    auto const lambda = _exponent;
    auto const a = _parameter_a;
    auto const b = _parameter_b;

    auto const s_eff = (s_L_within_range - s_L_res) / (s_L_max - s_L_res);

    auto const dk_rel_LRdsL = a * b * std::pow(1. - s_L_within_range, b - 1.);

    auto const _2L_L = (2. + lambda) / lambda;
    auto const s_G_eff = 1. - s_eff;
    auto const dk_rel_GRdse =
        -2. * s_G_eff * (1. - std::pow(s_eff, _2L_L)) -
        _2L_L * std::pow(s_eff, _2L_L - 1.) * s_G_eff * s_G_eff;
    auto const dk_rel_GRdsL = dk_rel_GRdse * _dse_dsL;

    return Eigen::Vector2d{dk_rel_LRdsL, dk_rel_GRdsL};
}

}  // namespace MaterialPropertyLib
