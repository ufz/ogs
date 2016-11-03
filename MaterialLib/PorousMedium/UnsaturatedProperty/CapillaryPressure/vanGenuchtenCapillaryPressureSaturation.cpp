/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   vanGenuchtenCapillaryPressureSaturation.cpp
 *
 *  Created on October 28, 2016, 6:05 PM
 */

#include "vanGenuchtenCapillaryPressureSaturation.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double vanGenuchtenCapillaryPressureSaturation::getCapillaryPressure(
    const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    const double pc =
        _pb * std::pow(std::pow(Se, (-1.0 / _mm)) - 1.0, 1.0 - _mm);
    return MathLib::limitValueInInterval(pc, _minor_offset, _Pc_max);
}

double vanGenuchtenCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc = (capillary_pressure < 0.0) ? 0.0 : capillary_pressure;
    double Se = std::pow(pc / _pb, 1.0 / (1.0 - _mm)) + 1.0;
    Se = std::pow(Se, -_mm);
    const double S = Se * (_Smax - _Sr) + _Sr;
    return MathLib::limitValueInInterval(S, _Sr + _minor_offset,
                                         _Smax - _minor_offset);
}

double vanGenuchtenCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double val1 = std::pow(((S - _Sr) / (_Smax - _Sr)), -1.0 / _mm);
    const double val2 = std::pow(val1 - 1.0, -_mm);
    return _pb * (_mm - 1.0) * val1 * val2 / (_mm * (S - _Sr));
}

}  // end namespace
}  // end namespace
