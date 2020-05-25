/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
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
    const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double pc = pb_ * std::pow(Se, -1.0 / m_);
    return std::clamp(pc, minor_offset_, pc_max_);
}

double BrooksCoreyCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc =
        (capillary_pressure < 0.0) ? minor_offset_ : capillary_pressure;
    const double Se = std::pow(pc / pb_, -m_);
    const double S = Se * (saturation_max_ - saturation_r_) + saturation_r_;
    return std::clamp(S, saturation_r_ + minor_offset_,
                      saturation_max_ - minor_offset_);
}

double BrooksCoreyCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double val = std::pow(
        ((S - saturation_r_) / (saturation_max_ - saturation_r_)), -1.0 / m_);
    return (pb_ * val) / (m_ * (saturation_r_ - S));
}

}  // namespace PorousMedium
}  // namespace MaterialLib
