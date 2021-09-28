/**
 * \file
 * \author Norbert Grunwald
 * \date   01.12.2020
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RelPermUdellNonwettingPhase.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermUdellNonwettingPhase::RelPermUdellNonwettingPhase(
    std::string name,
    const double residual_liquid_saturation,
    const double residual_gas_saturation,
    const double min_relative_permeability)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      min_relative_permeability_(min_relative_permeability)
{
    name_ = std::move(name);
}

PropertyDataType RelPermUdellNonwettingPhase::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (std::isnan(S_L))
    {
        OGS_FATAL(
            "In RelPermUdellNonwettingPhase::value, the liquid saturation is "
            "NaN.");
    }

    auto const S_L_res = residual_liquid_saturation_;
    auto const S_L_max = 1. - residual_gas_saturation_;

    auto const S_e = (S_L - S_L_res) / (S_L_max - S_L_res);

    if (S_e >= 1.0)
    {
        // fully saturated medium
        return min_relative_permeability_;
    }
    if (S_e <= 0.0)
    {
        // dry medium
        return 1.0;
    }

    auto const S_e_g = (1. - S_e);

    return std::max(min_relative_permeability_, S_e_g * S_e_g * S_e_g);
}
PropertyDataType RelPermUdellNonwettingPhase::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (primary_variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "RelPermUdellNonwettingPhase::dValue is implemented for "
            "derivatives with respect to liquid saturation only.");
    }

    const double S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    auto const S_L_res = residual_liquid_saturation_;
    auto const S_L_max = 1. - residual_gas_saturation_;
    auto const S_e = (S_L - S_L_res) / (S_L_max - S_L_res);

    if ((S_e < 0.) || (S_e > 1.))
    {
        return 0.;
    }

    auto const dS_e_dS_L = 1. / (S_L_max - S_L_res);

    auto const dk_rel_GR_dS_e = -3. * (1. - S_e) * (1. - S_e);
    return dk_rel_GR_dS_e * dS_e_dS_L;
}

}  // namespace MaterialPropertyLib
