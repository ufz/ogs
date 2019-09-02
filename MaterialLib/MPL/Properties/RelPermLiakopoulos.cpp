/**
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

#include "MaterialLib/MPL/Properties/RelPermLiakopoulos.h"
#include "MaterialLib/MPL/Medium.h"

#include <algorithm>
#include <cmath>

namespace MaterialPropertyLib
{
/**
 */
PropertyDataType RelPermLiakopoulos::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos,
    double const t) const
{
    /// here, an extra computation of saturation is forced, guaranteeing a
    /// correct value. In order to speed up the computing time, saturation could
    /// be insertred into the primary variable array after it is computed in the
    /// FEM assembly.
    auto const s_L = _medium->property(PropertyType::saturation)
                         .template value<double>(variable_array, pos, t);

    auto const s_L_res = _residual_liquid_saturation;
    auto const k_rel_min_GR = _min_relative_permeability_gas;
    auto const a = _parameter_a;
    auto const b = _parameter_b;
    auto const lambda = _exponent;

    auto const s_eff = (s_L - s_L_res) / (1. - s_L_res);

    auto const k_rel_LR = 1. - a * std::pow(1. - s_L, b);
    auto const k_rel_GR = (1. - s_eff) * (1. - s_eff) *
                          (1. - std::pow(s_eff, (2. + lambda) / lambda));

    const Pair kRel = {std::max(k_rel_LR, 0.),
                       std::max(k_rel_GR, k_rel_min_GR)};

    return kRel;
}
PropertyDataType RelPermLiakopoulos::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t) const
{
    assert((primary_variable == Variable::liquid_saturation) &&
           "RelPermLiakopoulos::dValue is implemented for "
           " derivatives with respect to liquid saturation only.");
    auto const s_L = _medium->property(PropertyType::saturation)
                         .template value<double>(variable_array, pos, t);

    auto const s_L_res = _residual_liquid_saturation;
    auto const s_L_max = 1.;
    auto const lambda = _exponent;
    auto const a = _parameter_a;
    auto const b = _parameter_b;

    auto const s_eff = (s_L - s_L_res) / (s_L_max - s_L_res);
    auto const d_se_d_sL = 1. / (s_L_max - s_L_res);

    auto const dk_rel_LRdsL = a * b * std::pow(1. - s_L, b - 1.);

    auto const _2L_L = (2. + lambda) / lambda;
    auto const dk_rel_GRdse =
        -2. * (1 - s_eff) * (1. - std::pow(s_eff, _2L_L)) -
        _2L_L * std::pow(s_eff, _2L_L - 1.) * (1. - s_eff) * (1. - s_eff);
    auto const dk_rel_GRdsL = dk_rel_GRdse * d_se_d_sL;

    const Pair dkReldsL = {{dk_rel_LRdsL, dk_rel_GRdsL}};

    return dkReldsL;
}

}  // namespace MaterialPropertyLib
