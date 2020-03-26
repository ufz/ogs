/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 26, 2020, 4:15 PM
 */

#include "CapillaryPressureLiakopoulos.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
PropertyDataType CapillaryPressureLiakopoulos::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S = std::clamp(saturation, _residual_liquid_saturation, 1.0);
    const double pc = std::pow((1.0 - S) / _parameter_a, 1.0 / _parameter_b);
    return std::clamp(pc, 0.0, _p_cap_max);
}

PropertyDataType CapillaryPressureLiakopoulos::dValue(
    VariableArray const& variable_array, Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S = std::clamp(saturation, _residual_liquid_saturation, 1.0);
    return -std::pow((1.0 - S) / _parameter_a, (1.0 / _parameter_b) - 1.0) /
           (_parameter_a * _parameter_b);
}

PropertyDataType CapillaryPressureLiakopoulos::d2Value(
    VariableArray const& variable_array, Variable const /*primary_variable1*/,
    Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S = std::clamp(saturation, _residual_liquid_saturation, 1.0);

    return ((1.0 / _parameter_b) - 1.0) *
           std::pow((1.0 - S) / _parameter_a, (1.0 / _parameter_b) - 2.0) /
           (_parameter_a * _parameter_a * _parameter_b);
}

}  // namespace MaterialPropertyLib
