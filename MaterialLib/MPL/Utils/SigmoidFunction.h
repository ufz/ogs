/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 20, 2022
 */

#pragma once

#include "MaterialLib/MPL/VariableType.h"  // for VariableArray
#include "ParameterLib/SpatialPosition.h"

namespace MaterialPropertyLib
{
/**
 * Sigmoid function with scalar argument and two constant parameters
 *
 * \details The sigmoid function is a smooth step function
 * with the range \f$[0,1]\f$ defined by the following formula
 *
 * \f[
 *      \left[1 + \exp(k(T - T_\mathrm{c})) \right]^{-1}
 * \f]
 *
 * where \f$T_\mathrm{c}\f$ is the critical value (e.g. a phase change
 * temperature). The parameter \f$k\f$ is proportional to the slope at
 * the characteristic value and controls thus the steepness. Letting
 * \f$k\f$ go to infinity, the heaviside step function is obtained.
 *
 */
class SigmoidFunction final
{
public:
    SigmoidFunction(double const k, double const T_c);

    double value(double const& T) const;

    double dValue(double const& T) const;

    double d2Value(double const& T) const;

private:
    double const k_;    //< steepness (slope parameter)
    double const T_c_;  //< characteristic value
                        // (location of the step)
};
}  // namespace MaterialPropertyLib
