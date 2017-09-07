/**
 * \file
 * \author Norbert Grunwald
 * \date   July 3rd, 2018
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pUndefined.h"
#include <string>
#include "MaterialLib/MPL/mpEnums.h"

namespace MaterialPropertyLib
{
Undefined::Undefined(PropertyEnum const& pEnum)
{
    thisPropertyEnum = pEnum;
}

PropertyDataType Undefined::value() const
{
    std::string property = convertEnumToString[thisPropertyEnum];

    OGS_FATAL(
        "The property \'%s\' (property-enum no. %i in "
        "MaterialLib/MPL/mpEnums.h) was requested, but is not defined in the "
        "project file.",
        property.c_str(), (int)thisPropertyEnum);
}

PropertyDataType Undefined::value(VariableArray const& /*v*/)
{
    return value();
}
}  // MaterialPropertyLib


