/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 *  Created on November 1, 2016, 9:45 AM
 */

#include "BrooksCoreyCapillaryPressureSaturation.h"

#include <algorithm>
#include <cmath>

namespace MaterialLib
{
namespace PorousMedium
{
double BrooksCoreyCapillaryPressureSaturation::getCapillaryPressure(
    const double saturation) const
{
    const double S = std::clamp(saturation, _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double pc = _pb * std::pow(Se, -1.0 / _m);
    return std::clamp(pc, _minor_offset, _pc_max);
}

double BrooksCoreyCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc =
        (capillary_pressure < 0.0) ? _minor_offset : capillary_pressure;
    const double Se = std::pow(pc / _pb, -_m);
    const double S = Se * (_saturation_max - _saturation_r) + _saturation_r;
    return std::clamp(S, _saturation_r + _minor_offset,
                      _saturation_max - _minor_offset);
}

double BrooksCoreyCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    const double S = std::clamp(saturation, _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double val = std::pow(
        ((S - _saturation_r) / (_saturation_max - _saturation_r)), -1.0 / _m);
    return (_pb * val) / (_m * (_saturation_r - S));
}

}  // namespace PorousMedium
}  // namespace MaterialLib
