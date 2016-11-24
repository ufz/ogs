/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 *  Created on November 1, 2016, 9:45 AM
 */

#include "BrookCoreyCapillaryPressureSaturation.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double BrookCoreyCapillaryPressureSaturation::getCapillaryPressure(
    const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    const double pc = _pb * std::pow(Se, -1.0 / _mm);
    return MathLib::limitValueInInterval(pc, _minor_offset, _Pc_max);
}

double BrookCoreyCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc =
        (capillary_pressure < 0.0) ? _minor_offset : capillary_pressure;
    const double Se = std::pow(pc / _pb, -_mm);
    const double S = Se * (_Smax - _Sr) + _Sr;
    return MathLib::limitValueInInterval(S, _Sr + _minor_offset,
                                         _Smax - _minor_offset);
}

double BrookCoreyCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double val = std::pow(((S - _Sr) / (_Smax - _Sr)), -1.0 / _mm);
    return (_pb * val) / (_mm * (_Sr - S));
}

}  // end namespace
}  // end namespace
