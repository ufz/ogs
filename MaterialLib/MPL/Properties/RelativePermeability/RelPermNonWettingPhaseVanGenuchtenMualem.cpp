/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 1, 2020, 2:15 PM
 */

#include "RelPermNonWettingPhaseVanGenuchtenMualem.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/CheckVanGenuchtenExponentRange.h"

namespace MaterialPropertyLib
{
RelPermNonWettingPhaseVanGenuchtenMualem::
    RelPermNonWettingPhaseVanGenuchtenMualem(const double Snr,
                                             const double Snmax, const double m,
                                             const double krel_min)
    : S_L_res_(1. - Snmax), S_L_max(1. - Snr), m_(m), krel_min_(krel_min)
{
    checkVanGenuchtenExponentRange(m_);
}

PropertyDataType RelPermNonWettingPhaseVanGenuchtenMualem::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        S_L_res_ + minor_offset_, S_L_max - minor_offset_);

    const double Se = (S_L - S_L_res_) / (S_L_max - S_L_res_);
    const double krel =
        std::sqrt(1.0 - Se) * std::pow(1.0 - std::pow(Se, 1.0 / m_), 2.0 * m_);
    return std::max(krel_min_, krel);
}

PropertyDataType RelPermNonWettingPhaseVanGenuchtenMualem::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert(
        (primary_variable == Variable::liquid_saturation) &&
        "RelPermNonWettingPhaseVanGenuchtenMualem::dValue is implemented for "
        "derivatives with respect to liquid saturation only.");

    const double S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    if (S_L < S_L_res_ + minor_offset_ || S_L > S_L_max - minor_offset_)
    {
        return 0.0;
    }

    const double Se = (S_L - S_L_res_) / (S_L_max - S_L_res_);
    const double val1 = std::sqrt(1.0 - Se);
    const double val2 = 1.0 - std::pow(Se, 1.0 / m_);
    return -0.5 * std::pow(val2, 2.0 * m_) / val1 -
           2.0 * std::pow(Se, -1.0 + 1.0 / m_) * val1 *
               std::pow(val2, 2.0 * m_ - 1.0);
}

}  // namespace MaterialPropertyLib