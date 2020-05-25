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

#include "NonWettingPhaseVanGenuchten.h"

#include <algorithm>
#include <cmath>

namespace MaterialLib
{
namespace PorousMedium
{
double NonWettingPhaseVanGenuchten::getValue(const double saturation_w) const
{
    const double S = std::clamp(saturation_w,
                                saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double krel = std::cbrt(1.0 - Se) *
                        std::pow(1.0 - std::pow(Se, 1.0 / m_), 2.0 * m_);
    return std::max(krel_min_, krel);
}

double NonWettingPhaseVanGenuchten::getdValue(const double saturation_w) const
{
    const double S = std::clamp(saturation_w,
                                saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double cbrt1_Se = std::cbrt(1.0 - Se);
    const double temp_val = 1.0 - std::pow(Se, 1.0 / m_);
    return (-std::pow(temp_val, 2. * m_) / (3. * cbrt1_Se * cbrt1_Se) -
            2. * cbrt1_Se * std::pow(temp_val, 2. * m_ - 1.) *
                std::pow(Se, (1. - m_) / m_)) /
           (saturation_max_ - saturation_r_);
}
}  // namespace PorousMedium
}  // namespace MaterialLib
