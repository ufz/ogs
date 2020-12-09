/**
 * \file
 * \author Norbert Grunwald
 * \date   01.12.2020
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RelPermUdell.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermUdell::RelPermUdell(std::string name,
                           const double residual_liquid_saturation,
                           const double residual_gas_saturation,
                           const double min_relative_permeability_liquid,
                           const double min_relative_permeability_gas)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      min_relative_permeability_liquid_(min_relative_permeability_liquid),
      min_relative_permeability_gas_(min_relative_permeability_gas)
{
    name_ = std::move(name);
}

PropertyDataType RelPermUdell::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double s_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (std::isnan(s_L))
    {
        OGS_FATAL("Liquid saturation not set in RelPermUdell::value().");
    }

    auto const s_L_res = residual_liquid_saturation_;
    auto const s_L_max = 1. - residual_gas_saturation_;
    auto const k_rel_min_LR = min_relative_permeability_liquid_;
    auto const k_rel_min_GR = min_relative_permeability_gas_;

    auto const s = (s_L - s_L_res) / (s_L_max - s_L_res);

    if (s >= 1.0)
    {
        // fully saturated medium
        return Eigen::Vector2d{1.0, k_rel_min_GR};
    }
    if (s <= 0.0)
    {
        // dry medium
        return Eigen::Vector2d{k_rel_min_LR, 1.0};
    }

    auto const k_rel_LR = s * s * s;
    auto const k_rel_GR = (1. - s) * (1. - s) * (1. - s);

    return Eigen::Vector2d{std::max(k_rel_LR, k_rel_min_LR),
                           std::max(k_rel_GR, k_rel_min_GR)};
}
PropertyDataType RelPermUdell::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "RelPermUdell::dValue is implemented for "
           " derivatives with respect to liquid saturation only.");

    const double s_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    auto const s_L_res = residual_liquid_saturation_;
    auto const s_L_max = 1. - residual_gas_saturation_;
    auto const s = (s_L - s_L_res) / (s_L_max - s_L_res);

    if ((s < 0.) || (s > 1.))
    {
        return Eigen::Vector2d{0., 0.};
    }

    auto const d_se_d_sL = 1. / (s_L_max - s_L_res);

    auto const dk_rel_LRdse = 3. * s * s;
    auto const dk_rel_LRdsL = dk_rel_LRdse * d_se_d_sL;

    auto const dk_rel_GRdse = -3. * (1. - s) * (1. - s);
    auto const dk_rel_GRdsL = dk_rel_GRdse * d_se_d_sL;

    return Eigen::Vector2d{dk_rel_LRdsL, dk_rel_GRdsL};
}

}  // namespace MaterialPropertyLib
