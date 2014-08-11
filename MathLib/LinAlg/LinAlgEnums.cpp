/**
 * \author Wenqing Wang
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinAlgEnums.h"

namespace MathLib
{

std::string convertVecNormTypeToString(VecNormType normType)
{
    switch (normType)
    {
    case VecNormType::NORM1: return "NORM1";
    case VecNormType::NORM2: return "NORM2";
    case VecNormType::INFINITY_N: return "INFINITY_N";
    default: return "INVALID";
    }
}

VecNormType convertVecNormTypeToString(const std::string &str)
{
    if (str.compare("NORM1")==0) return VecNormType::NORM1;
    if (str.compare("NORM2")==0) return VecNormType::NORM2;
    if (str.compare("INFINITY_N")==0) return VecNormType::INFINITY_N;
    return VecNormType::INVALID;
}

} // end namespace MathLib
