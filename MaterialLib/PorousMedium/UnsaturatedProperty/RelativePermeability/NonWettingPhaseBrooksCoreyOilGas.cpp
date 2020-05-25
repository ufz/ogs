/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 2, 2016, 10:47 AM
 */

#include "NonWettingPhaseBrooksCoreyOilGas.h"

#include <algorithm>
#include <cmath>

namespace MaterialLib
{
namespace PorousMedium
{
double NonWettingPhaseBrooksCoreyOilGas::getValue(
    const double saturation_w) const
{
    const double S = std::clamp(saturation_w,
                                saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double krel =
        (1.0 - Se) * (1.0 - Se) * (1.0 - std::pow(Se, 1.0 + 2.0 / m_));
    return std::max(krel_min_, krel);
}

double NonWettingPhaseBrooksCoreyOilGas::getdValue(
    const double saturation_w) const
{
    const double S = std::clamp(saturation_w,
                                saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    return (-2. * (1.0 - Se) * (1.0 - std::pow(Se, 1.0 + 2.0 / m_)) -
            (1.0 + 2.0 / m_) * (1.0 - Se) * (1.0 - Se) *
                std::pow(Se, 2.0 / m_)) /
           (saturation_max_ - saturation_r_);
}

}  // namespace PorousMedium
}  // namespace MaterialLib
