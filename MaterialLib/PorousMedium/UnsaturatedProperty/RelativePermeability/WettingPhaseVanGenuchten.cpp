/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WettingPhaseVanGenuchten.cpp
 *
 * Created on November 2, 2016, 11:24 AM
 */

#include "WettingPhaseVanGenuchten.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double WettingPhaseVanGenuchten::getValue(const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    const double val = 1.0 - std::pow(1.0 - std::pow(Se, 1.0 / _mm), _mm);
    const double krel = std::sqrt(Se) * val * val;
    return std::max(_krel_min,  krel);
}

double WettingPhaseVanGenuchten::getdValue(const double saturation) const
{
    const double S = MathLib::limitValueInInterval(
        saturation, _Sr + _minor_offset, _Smax - _minor_offset);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    const double sqrtSe = std::sqrt(Se);
    const double temp_val = 1.0 - std::pow(1.0 - std::pow(Se, 1.0 / _mm), _mm);
    return (0.5 * temp_val * temp_val / sqrtSe +
            2. * sqrtSe * temp_val *
                std::pow(1.0 - std::pow(Se, 1.0 / _mm), _mm - 1.) *
                std::pow(Se, (1.0 - _mm) / _mm)) /
           (_Smax - _Sr);
}

}  // end namespace
}  // end namespace
