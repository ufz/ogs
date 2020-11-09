/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaturationExponential.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/MathTools.h"

namespace MaterialPropertyLib
{
SaturationExponential::SaturationExponential(
    std::string name,
    const double residual_liquid_saturation,
    const double residual_gas_saturation,
    const double p_cap_ref,
    const double exponent)
    : residual_liquid_saturation_(residual_liquid_saturation),
      residual_gas_saturation_(residual_gas_saturation),
      p_cap_ref_(p_cap_ref),
      exponent_(exponent)
{
    name_ = std::move(name);
    if (p_cap_ref_ <= 0.)
    {
        OGS_FATAL(
            "Reference capillary pressure must be positive in "
            "MPL::SaturationExponential.");
    }
};

PropertyDataType SaturationExponential::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    const double s_res = residual_liquid_saturation_;
    const double s_max = 1. - residual_gas_saturation_;

    const double pc = std::clamp(p_cap, 0., p_cap_ref_);
    const double s_e = 1. - std::pow(pc / p_cap_ref_, exponent_);
    return s_e * (s_max - s_res) + s_res;
}

PropertyDataType SaturationExponential::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::capillary_pressure) &&
           "SaturationExponential::dValue is implemented for derivatives with "
           "respect to capillary pressure only.");

    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    const double s_res = residual_liquid_saturation_;
    const double s_max = 1. - residual_gas_saturation_;

    if ((p_cap > p_cap_ref_) || (p_cap <= 0.))
    {
        return 0.;
    }
    return (exponent_ / p_cap) * (s_res - s_max) *
           std::pow(p_cap / p_cap_ref_, exponent_);
}

PropertyDataType SaturationExponential::d2Value(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable1*/, Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("SaturationExponential::d2Value() is not implemented.");
}

}  // namespace MaterialPropertyLib
