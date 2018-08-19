/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaturationVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
SaturationVanGenuchten::SaturationVanGenuchten(
    double const residual_liquid_saturation,
    double const residual_gas_saturation,
    double const exponent,
    double const entry_pressure)
    : _S_L_res(residual_liquid_saturation),
      _S_L_max(1. - residual_gas_saturation),
      _m(exponent),
      _p_b(entry_pressure)
{
    if (!(_m > 0 && _m < 1))
    {
        OGS_FATAL(
            "The exponent value m = %g of van Genuchten saturation model, is "
            "out of its range of (0, 1)",
            _m);
    }
}

PropertyDataType SaturationVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/) const
{
    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if (p_cap <= 0)
    {
        return _S_L_max;
    }

    double const p = p_cap / _p_b;
    double const n = 1. / (1. - _m);
    double const p_to_n = std::pow(p, n);

    double const S_eff = std::pow(p_to_n + 1., -_m);
    double const S = S_eff * _S_L_max - S_eff * _S_L_res + _S_L_res;
    return std::clamp(S, _S_L_res, _S_L_max);
}

PropertyDataType SaturationVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::capillary_pressure) &&
           "SaturationVanGenuchten::dvalue is implemented for derivatives with "
           "respect to capillary pressure only.");

    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if (p_cap <= 0)
    {
        return 0;
    }

    double const p = p_cap / _p_b;
    double const n = 1. / (1. - _m);
    double const p_to_n = std::pow(p, n);

    double const S_eff = std::pow(p_to_n + 1., -_m);
    double const S = S_eff * _S_L_max - S_eff * _S_L_res + _S_L_res;

    if (S < _S_L_res || S > _S_L_max)
    {
        return 0;
    }

    double const dS_eff_dp_cap = -_m * std::pow(p_cap / _p_b, n - 1) *
                                 std::pow(1 + p_to_n, -1. - _m) /
                                 (_p_b * (1. - _m));
    return dS_eff_dp_cap * (_S_L_max - _S_L_res);
}

PropertyDataType SaturationVanGenuchten::d2Value(
    VariableArray const& variable_array, Variable const primary_variable1,
    Variable const primary_variable2,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/) const
{
    (void)primary_variable1;
    (void)primary_variable2;
    assert((primary_variable1 == Variable::capillary_pressure) &&
           (primary_variable2 == Variable::capillary_pressure) &&
           "SaturationVanGenuchten::d2Value is implemented for "
           " derivatives with respect to capillary pressure only.");

    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if (p_cap <= 0)
    {
        return 0;
    }

    double const p = p_cap / _p_b;
    double const n = 1. / (1. - _m);
    double const p_to_n = std::pow(p, n);

    double const S_eff = std::pow(p_to_n + 1., -_m);
    double const S = S_eff * _S_L_max - S_eff * _S_L_res + _S_L_res;

    if (S < _S_L_res || S > _S_L_max)
    {
        return 0;
    }

    double const d2S_eff_dp_cap2 =
        _m * p_to_n * std::pow(p_to_n + 1., -_m - 2.) * (p_to_n - _m) /
        (p_cap * p_cap * (_m - 1.) * (_m - 1.));
    return d2S_eff_dp_cap2 * (_S_L_max - _S_L_res);
}
}  // namespace MaterialPropertyLib
