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
        saturation, _Sr + _perturbation, _Smax - _perturbation);
    const double Se = (S - _Sr) / (_Smax - _Sr);
    const double Krel =
        std::sqrt(Se) *
        std::pow(1.0 - std::pow(1.0 - std::pow(Se, 1.0 / _mm), _mm), 2);
    return Krel < _Krel_min ? _Krel_min : Krel;
}

}  // end namespace
}  // end namespace
