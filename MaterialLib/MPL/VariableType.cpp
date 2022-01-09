/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VariableType.h"

#include <boost/algorithm/string/predicate.hpp>

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
Variable convertStringToVariable(std::string const& string)
{
    for (int i = 0; i < static_cast<int>(Variable::number_of_variables); ++i)
    {
        if (boost::iequals(string, variable_enum_to_string[i]))
        {
            return static_cast<Variable>(i);
        }
    }

    OGS_FATAL(
        "The variable name '{:s}' does not correspond to any known variable",
        string);
}
}  // namespace MaterialPropertyLib
