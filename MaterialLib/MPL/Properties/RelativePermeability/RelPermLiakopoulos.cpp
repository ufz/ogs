/**
 * \file
 * \author Norbert Grunwald
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
RelPermLiakopoulos::RelPermLiakopoulos(std::string name)
{
    name_ = std::move(name);
}

PropertyDataType RelPermLiakopoulos::value(
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

    if (s_L <= s_L_res)
    {
        return 0.0;
    }

    if (s_L >= 1.)
    {
        return 1.0;
    }

    auto const a = parameter_a_;
    auto const b = parameter_b_;

    auto const k_rel_LR = 1. - a * std::pow(1. - s_L, b);

    return std::max(k_rel_LR, 0.);
}

PropertyDataType RelPermLiakopoulos::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "RelPermLiakopoulos::dValue is implemented for derivatives with "
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
    auto const s_L_max = maximal_liquid_saturation_;

    const double s_L_within_range = std::min(std::max(s_L_res, s_L), s_L_max);

    auto const a = parameter_a_;
    auto const b = parameter_b_;

    return a * b * std::pow(1. - s_L_within_range, b - 1.);
}

}  // namespace MaterialPropertyLib
