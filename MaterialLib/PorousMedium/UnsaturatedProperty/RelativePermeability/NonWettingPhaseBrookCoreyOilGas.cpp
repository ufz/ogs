/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   NonWettingPhaseBrookCoreyOilGas.cpp
 *
 * Created on November 2, 2016, 10:47 AM
 */

#include "NonWettingPhaseBrookCoreyOilGas.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double NonWettingPhaseBrookCoreyOilGas::getValue(
    const double saturation_w) const
{
    const double S = MathLib::limitValueInInterval(
        saturation_w, _Sr + _minor_offset, _Smax - _minor_offset);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    const double krel =
        (1.0 - Se) * (1.0 - Se) * (1.0 - std::pow(Se, 1.0 + 2.0 / _mm));
    return std::max(_krel_min,  krel);
}

double NonWettingPhaseBrookCoreyOilGas::getdValue(
    const double saturation_w) const
{
    const double S = MathLib::limitValueInInterval(
        saturation_w, _Sr + _minor_offset, _Smax - _minor_offset);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    return (-2. * (1.0 - Se) * (1.0 - std::pow(Se, 1.0 + 2.0 / _mm)) -
            (1.0 + 2.0 / _mm) * (1.0 - Se) * (1.0 - Se) *
                std::pow(Se, 2.0 / _mm)) /
           (_Smax - _Sr);
}

}  // end namespace
}  // end namespace
