/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
                                saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double val = 1.0 - std::pow(1.0 - std::pow(Se, 1.0 / m_), m_);
    const double krel = std::sqrt(Se) * val * val;
    return std::max(krel_min_, krel);
}

double WettingPhaseVanGenuchten::getdValue(const double saturation) const
{
    const double S = std::clamp(saturation,
                                saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double sqrtSe = std::sqrt(Se);
    const double temp_val = 1.0 - std::pow(1.0 - std::pow(Se, 1.0 / m_), m_);
    return (0.5 * temp_val * temp_val / sqrtSe +
            2. * sqrtSe * temp_val *
                std::pow(1.0 - std::pow(Se, 1.0 / m_), m_ - 1.) *
                std::pow(Se, (1.0 - m_) / m_)) /
           (saturation_max_ - saturation_r_);
}

}  // namespace PorousMedium
}  // namespace MaterialLib
