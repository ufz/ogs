/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 20, 2020, 9:59 AM
 */

#include "CapillaryPressureVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
CapillaryPressureVanGenuchten::CapillaryPressureVanGenuchten(
    double const residual_liquid_saturation,
    double const maximum_liquid_saturation, double const exponent,
    double const p_b, double const max_capillary_pressure)
    : _residual_saturation(residual_liquid_saturation),
      _maximuml_saturation(maximum_liquid_saturation),
      _m(exponent),
      _p_b(p_b),
      _pc_max(max_capillary_pressure)
{
    if (!(_m > 0 && _m < 1))
    {
        OGS_FATAL(
            "The exponent value m = %g of van Genuchten saturation model, is "
            "out of its range of (0, 1)",
            _m);
    }
}

PropertyDataType CapillaryPressureVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);
    const double Se = (S - _residual_saturation) /
                      (_maximuml_saturation - _residual_saturation);
    const double pc =
        (Se < 1.0) ? _p_b * std::pow(std::pow(Se, (-1.0 / _m)) - 1.0, 1.0 - _m)
                   : 0.0;
    return std::clamp(pc, 0.0, _pc_max);
}

PropertyDataType CapillaryPressureVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);
    const double Se = (S - _residual_saturation) /
                      (_maximuml_saturation - _residual_saturation);

    if (!(Se < 1.0))
    {
        return 0.0;
    }
    const double val1 = std::pow(Se, -1.0 / _m);
    const double val2 = std::pow(val1 - 1.0, -_m);
    return _p_b * (_m - 1.0) * val1 * val2 / (_m * (S - _residual_saturation));
}

PropertyDataType CapillaryPressureVanGenuchten::d2Value(
    VariableArray const& variable_array, Variable const /*primary_variable1*/,
    Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double saturation = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S =
        std::clamp(saturation, _residual_saturation, _maximuml_saturation);
    const double Se = (S - _residual_saturation) /
                      (_maximuml_saturation - _residual_saturation);
    if (!(Se < 1.0))
    {
        return 0.0;
    }

    const double val1 = std::pow(Se, 1.0 / _m);
    return -_p_b /
           (_m * _m * (S - _residual_saturation) * (S - _residual_saturation)) *
           std::pow(1 - val1, -_m - 1) * std::pow(val1, _m - 1) *
           ((1 - _m * _m) * val1 + _m - 1);
}

}  // namespace MaterialPropertyLib
