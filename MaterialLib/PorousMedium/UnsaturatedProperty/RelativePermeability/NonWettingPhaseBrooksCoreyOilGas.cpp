/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 2, 2016, 10:47 AM
 */

#include "NonWettingPhaseBrooksCoreyOilGas.h"

#include <algorithm>
#include <cmath>

namespace MaterialLib
{
namespace PorousMedium
{
double NonWettingPhaseBrooksCoreyOilGas::getValue(
    const double saturation_w) const
{
    const double S = std::clamp(saturation_w,
                                _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double krel =
        (1.0 - Se) * (1.0 - Se) * (1.0 - std::pow(Se, 1.0 + 2.0 / _m));
    return std::max(_krel_min, krel);
}

double NonWettingPhaseBrooksCoreyOilGas::getdValue(
    const double saturation_w) const
{
    const double S = std::clamp(saturation_w,
                                _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    return (-2. * (1.0 - Se) * (1.0 - std::pow(Se, 1.0 + 2.0 / _m)) -
            (1.0 + 2.0 / _m) * (1.0 - Se) * (1.0 - Se) *
                std::pow(Se, 2.0 / _m)) /
           (_saturation_max - _saturation_r);
}

}  // namespace PorousMedium
}  // namespace MaterialLib
