/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RelPermVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermVanGenuchten::RelPermVanGenuchten(
    double const residual_liquid_saturation,
    double const residual_gas_saturation,
    double const min_relative_permeability_liquid,
    double const exponent)
    : _S_L_res(residual_liquid_saturation),
      _S_L_max(1. - residual_gas_saturation),
      _k_rel_min(min_relative_permeability_liquid),
      _m(exponent)
{
    if (!(_m > 0 && _m < 1))
    {
        OGS_FATAL(
            "The exponent value m = %g of van Genuchten relative permeability "
            "model, is out of its range of (0, 1)",
            _m);
    }
}

PropertyDataType RelPermVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/) const
{
    double const S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        _S_L_res, _S_L_max);

    double const S_eff = (S_L - _S_L_res) / (_S_L_max - _S_L_res);
    double const v = 1. - std::pow(1. - std::pow(S_eff, 1. / _m), _m);
    double const k_rel = std::sqrt(S_eff) * v * v;
    return std::max(_k_rel_min, k_rel);
}

PropertyDataType RelPermVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "RelativePermeabilityVanGenuchten::dValue is implemented for "
           "derivatives with respect to liquid saturation only.");
    double const S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        _S_L_res, _S_L_max);

    double const S_eff = (S_L - _S_L_res) / (_S_L_max - _S_L_res);
    if (S_eff <= 0)  // prevent division by zero
    {
        return 0;
    }

    if (S_eff >= 1)  // prevent taking root of zero
    {
        return 0;
    }

    double const S_eff_to_1_over_m = std::pow(S_eff, 1. / _m);
    double const v = 1. - std::pow(1. - S_eff_to_1_over_m, _m);
    double const sqrt_S_eff = std::sqrt(S_eff);
    double const k_rel = sqrt_S_eff * v * v;

    if (k_rel < _k_rel_min)
    {
        return 0;
    }

    return (0.5 * v * v / sqrt_S_eff +
            2. * sqrt_S_eff * v * std::pow(1. - S_eff_to_1_over_m, _m - 1.) *
                S_eff_to_1_over_m / S_eff) /
           (_S_L_max - _S_L_res);
}

}  // namespace MaterialPropertyLib
