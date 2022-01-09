/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 2, 2016, 11:24 AM
 */

#include "WettingPhaseVanGenuchten.h"

#include <algorithm>
#include <cmath>

namespace MaterialLib
{
namespace PorousMedium
{
double WettingPhaseVanGenuchten::getValue(const double saturation) const
{
    const double S = std::clamp(saturation,
                                _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double val = 1.0 - std::pow(1.0 - std::pow(Se, 1.0 / _m), _m);
    const double krel = std::sqrt(Se) * val * val;
    return std::max(_krel_min, krel);
}

double WettingPhaseVanGenuchten::getdValue(const double saturation) const
{
    const double S = std::clamp(saturation,
                                _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double sqrtSe = std::sqrt(Se);
    const double temp_val = 1.0 - std::pow(1.0 - std::pow(Se, 1.0 / _m), _m);
    return (0.5 * temp_val * temp_val / sqrtSe +
            2. * sqrtSe * temp_val *
                std::pow(1.0 - std::pow(Se, 1.0 / _m), _m - 1.) *
                std::pow(Se, (1.0 - _m) / _m)) /
           (_saturation_max - _saturation_r);
}

}  // namespace PorousMedium
}  // namespace MaterialLib
