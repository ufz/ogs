/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaturationBrooksCorey.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/MathTools.h"

namespace MaterialPropertyLib
{
SaturationBrooksCorey::SaturationBrooksCorey(
    std::string name,
    const double residual_liquid_saturation,
    const double residual_gas_saturation,
    const double exponent,
    const double entry_pressure)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      exponent_(exponent),
      entry_pressure_(entry_pressure)
{
    name_ = std::move(name);
};

PropertyDataType SaturationBrooksCorey::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double p_cap = variable_array.capillary_pressure;

    const double s_L_res = residual_liquid_saturation_;
    const double s_L_max = 1.0 - residual_gas_saturation_;
    const double lambda = exponent_;
    const double p_b = entry_pressure_;

    if (p_cap <= p_b)
        return s_L_max;

    const double s_eff = std::pow(p_b / p_cap, lambda);
    return s_eff * (s_L_max - s_L_res) + s_L_res;
}

PropertyDataType SaturationBrooksCorey::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "SaturationBrooksCorey::dValue is implemented for derivatives with "
            "respect to capillary pressure only.");
    }

    const double p_b = entry_pressure_;
    const double p_cap = variable_array.capillary_pressure;

    if (p_cap <= p_b)
    {
        return 0.;
    }

    const double s_L_res = residual_liquid_saturation_;
    const double s_L_max = 1.0 - residual_gas_saturation_;
    const double lambda = exponent_;

    const double ds_eff_dp_cap =
        -lambda * std::pow(p_b, lambda) / std::pow(p_cap, lambda + 1);
    return ds_eff_dp_cap * (s_L_max - s_L_res);
}

PropertyDataType SaturationBrooksCorey::d2Value(
    VariableArray const& variable_array, Variable const variable1,
    Variable const variable2, ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/, double const /*dt*/) const
{
    if ((variable1 != Variable::capillary_pressure) &&
        (variable2 != Variable::capillary_pressure))
    {
        OGS_FATAL(
            "SaturationBrooksCorey::d2Value is implemented for derivatives "
            "with respect to capillary pressure only.");
    }

    const double p_b = entry_pressure_;
    const double p_cap = std::max(p_b, variable_array.capillary_pressure);

    if (p_cap <= p_b)
    {
        return 0.;
    }

    const double s_L_res = residual_liquid_saturation_;
    const double s_L_max = 1.0 - residual_gas_saturation_;

    const double lambda = exponent_;

    return lambda * (lambda + 1) * std::pow(p_b / p_cap, lambda) /
           (p_cap * p_cap) * (s_L_max - s_L_res);
}

}  // namespace MaterialPropertyLib
