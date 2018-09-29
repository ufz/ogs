/**
 * \author Wenqing Wang
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

namespace MathLib
{

/// Norm type. Not declared as class type in order to use the members as integers.
enum class VecNormType
{
    NORM1,        ///< \f$\sum_i |x_i|\f$
    NORM2,        ///< \f$\sqrt(\sum_i (x_i)^2)\f$
    INFINITY_N,    ///< \f$\mathrm{max}_i |x_i|\f$
    INVALID
};

/// convert VecNormType to string
std::string convertVecNormTypeToString(VecNormType normType);

/// convert string to VecNormType
VecNormType convertStringToVecNormType(const std::string &str);

} // end namespace MathLib
