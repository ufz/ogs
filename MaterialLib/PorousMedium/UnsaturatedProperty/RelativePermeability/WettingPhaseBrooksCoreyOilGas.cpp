/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WettingPhaseBrooksCoreyOilGas.cpp
 *
 * Created on November 1, 2016, 3:37 PM
 */

#include "WettingPhaseBrooksCoreyOilGas.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
double WettingPhaseBrooksCoreyOilGas::getValue(const double saturation) const
{
    const double S =
        MathLib::limitValueInInterval(saturation, _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double krel = std::pow(Se, 3.0 + 2.0 / _m);
    return std::max(_krel_min, krel);
}

double WettingPhaseBrooksCoreyOilGas::getdValue(const double saturation) const
{
    const double S =
        MathLib::limitValueInInterval(saturation, _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    return ((3.0 + 2.0 / _m) * std::pow(Se, 2.0 + 2.0 / _m)) /
           (_saturation_max - _saturation_r);
}
}  // end namespace
}  // end namespace
