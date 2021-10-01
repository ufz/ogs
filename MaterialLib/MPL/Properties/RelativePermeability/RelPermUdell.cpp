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

#include "RelPermUdell.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermUdell::RelPermUdell(std::string name,
                           const double residual_liquid_saturation,
                           const double residual_gas_saturation,
                           const double min_relative_permeability)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      min_relative_permeability_(min_relative_permeability)
{
    name_ = std::move(name);
}

PropertyDataType RelPermUdell::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (std::isnan(S_L))
    {
        OGS_FATAL("Liquid saturation not set in RelPermUdell::value().");
    }

    auto const S_L_res = residual_liquid_saturation_;
    auto const S_L_max = 1. - residual_gas_saturation_;

    auto const S_e = (S_L - S_L_res) / (S_L_max - S_L_res);

    if (S_e >= 1.0)
    {
        // fully saturated medium
        return 1.0;
    }
    if (S_e <= 0.0)
    {
        // dry medium
        return min_relative_permeability_;
    }

    return std::max(min_relative_permeability_, S_e * S_e * S_e);
}
PropertyDataType RelPermUdell::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "RelPermUdell::dValue is implemented for derivatives with respect "
            "to liquid saturation only.");
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

    auto const dk_rel_LR_dS_e = 3. * S_e * S_e;
    return dk_rel_LR_dS_e * dS_e_dS_L;
}

}  // namespace MaterialPropertyLib
