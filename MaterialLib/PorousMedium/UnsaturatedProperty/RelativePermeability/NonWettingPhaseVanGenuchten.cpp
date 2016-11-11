/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   NonWettingPhaseVanGenuchten.cpp
 *
 * Created on November 2, 2016, 11:24 AM
 */

#include "NonWettingPhaseVanGenuchten.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double NonWettingPhaseVanGenuchten::getValue(const double saturation_w) const
{
    const double S =
        MathLib::limitValueInInterval(saturation_w,
                                      _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double krel = std::cbrt(1.0 - Se) *
                        std::pow(1.0 - std::pow(Se, 1.0 / _mm), 2.0 * _mm);
    return std::max(_krel_min, krel);
}

double NonWettingPhaseVanGenuchten::getdValue(const double saturation_w) const
{
    const double S =
        MathLib::limitValueInInterval(saturation_w,
                                      _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double cbrt1_Se = std::cbrt(1.0 - Se);
    const double temp_val = 1.0 - std::pow(Se, 1.0 / _mm);
    return (-std::pow(temp_val, 2. * _mm) / (3. * cbrt1_Se * cbrt1_Se) -
            2. * cbrt1_Se * std::pow(temp_val, 2. * _mm - 1.) *
                std::pow(Se, (1. - _mm) / _mm)) /
           (_saturation_max - _saturation_r);
}
}  // end namespace
}  // end namespace
