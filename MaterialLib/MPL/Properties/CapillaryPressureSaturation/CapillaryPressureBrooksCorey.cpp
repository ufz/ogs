/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 25, 2020, 1:51 PM
 */

#include "CapillaryPressureBrooksCorey.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
PropertyDataType CapillaryPressureBrooksCorey::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S = std::clamp(saturation, _residual_liquid_saturation,
                                _max_liquid_saturation);
    const double Se = (S - _residual_liquid_saturation) /
                      (_max_liquid_saturation - _residual_liquid_saturation);
    const double pc = _entry_pressure * std::pow(Se, -1.0 / _exponent);
    return std::clamp(pc, 0.0, _pc_max);
}

PropertyDataType CapillaryPressureBrooksCorey::dValue(
    VariableArray const& variable_array, Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S = std::clamp(saturation, _residual_liquid_saturation,
                                _max_liquid_saturation);
    const double Se = (S - _residual_liquid_saturation) /
                      (_max_liquid_saturation - _residual_liquid_saturation);

    return -_entry_pressure *std::pow(Se, -1.0 / _exponent) /
           (_exponent * (S - _residual_liquid_saturation));
}

PropertyDataType CapillaryPressureBrooksCorey::d2Value(
    VariableArray const& variable_array, Variable const /*primary_variable1*/,
    Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S = std::clamp(saturation, _residual_liquid_saturation,
                                _max_liquid_saturation);

    const double fac = (S - _residual_liquid_saturation);
    const double Se =
        fac / (_max_liquid_saturation - _residual_liquid_saturation);

    return _entry_pressure * (_exponent + 1.) * std::pow(Se, -1.0 / _exponent) /
           (_exponent * _exponent * fac * fac);
}

}  // namespace MaterialPropertyLib
