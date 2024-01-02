/**
 * \file
 * \author Norbert Grunwald
 * \date   02.07.2018
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
    if ((s_eff < 0.) || (s_eff > 1.))
    {
        return 0.;
    }

    auto const d_se_d_sL = 1. / (s_L_max - s_L_res);
    auto const dk_rel_LRdse =
        (3 * lambda + 2.) / lambda * std::pow(s_eff, 2. / lambda + 2.);

    return dk_rel_LRdse * d_se_d_sL;
}

}  // namespace MaterialPropertyLib
