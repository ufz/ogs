/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 1, 2016, 3:37 PM
 */

#include "WettingPhaseBrooksCoreyOilGas.h"

#include <algorithm>
#include <cmath>

namespace MaterialLib
{
namespace PorousMedium
{
double WettingPhaseBrooksCoreyOilGas::getValue(const double saturation) const
{
    const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double krel = std::pow(Se, 3.0 + 2.0 / m_);
    return std::max(krel_min_, krel);
}

double WettingPhaseBrooksCoreyOilGas::getdValue(const double saturation) const
{
    const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    return ((3.0 + 2.0 / m_) * std::pow(Se, 2.0 + 2.0 / m_)) /
           (saturation_max_ - saturation_r_);
}
}  // namespace PorousMedium
}  // namespace MaterialLib
