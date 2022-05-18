/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 16, 2019, 3:40 PM
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
 * with the range [0,1] defined by the following formula
 * 
 * \f[
 *      \left[1 + \exp(k(T - T_\text{c})) \right]^{-1}
 * \f]
 * 
 * where \f$T_\text{c})\f$ is the critical value 
 * (e.g. a phase change temperature).
 */
class SigmoidFunction final
{
public:
    SigmoidFunction(double const k, double const T_c);

    double value(double const& T) const;

    double dValue(double const& T) const;

    double d2Value(double const& T) const;

private:
    double const k_;    //< factor controlling the spreading
    double const T_c_;  //< characteristic value 
                        // (location of the step)
};
}  // namespace MaterialPropertyLib
