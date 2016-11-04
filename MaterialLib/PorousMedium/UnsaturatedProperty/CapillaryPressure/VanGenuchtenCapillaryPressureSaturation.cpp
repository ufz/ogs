/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   VanGenuchtenCapillaryPressureSaturation.cpp
 *
 *  Created on October 28, 2016, 6:05 PM
 */

#include "VanGenuchtenCapillaryPressureSaturation.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double VanGenuchtenCapillaryPressureSaturation::getCapillaryPressure(
    const double saturation) const
{
    const double S =
        MathLib::limitValueInInterval(saturation, _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double pc = _pb * std::pow(std::pow(Se, (-1.0 / _m)) - 1.0, 1.0 - _m);
    return MathLib::limitValueInInterval(pc, _minor_offset, _pc_max);
}

double VanGenuchtenCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc =
        (capillary_pressure < 0.0) ? _minor_offset : capillary_pressure;
    double Se = std::pow(pc / _pb, 1.0 / (1.0 - _m)) + 1.0;
    Se = std::pow(Se, -_m);
    const double S = Se * (_saturation_max - _saturation_r) + _saturation_r;
    return MathLib::limitValueInInterval(S, _saturation_r + _minor_offset,
                                         _saturation_max - _minor_offset);
}

double VanGenuchtenCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    const double S =
        MathLib::limitValueInInterval(saturation, _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double val1 = std::pow(
        ((S - _saturation_r) / (_saturation_max - _saturation_r)), -1.0 / _m);
    const double val2 = std::pow(val1 - 1.0, -_m);
    return _pb * (_m - 1.0) * val1 * val2 / (_m * (S - _saturation_r));
}

}  // end namespace
}  // end namespace
