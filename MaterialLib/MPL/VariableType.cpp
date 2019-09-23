/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VariableType.h"
#include <boost/algorithm/string/predicate.hpp>
#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
Variable convertStringToVariable(std::string const& input)
{
    if (boost::iequals(input, "concentration"))
    {
        return Variable::concentration;
    }
    if (boost::iequals(input, "phase_pressure"))
    {
        return Variable::phase_pressure;
    }
    if (boost::iequals(input, "capillary_pressure"))
    {
        return Variable::capillary_pressure;
    }
    if (boost::iequals(input, "density"))
    {
        return Variable::density;
    }
    if (boost::iequals(input, "temperature"))
    {
        return Variable::temperature;
    }
    if (boost::iequals(input, "liquid_saturation"))
    {
        return Variable::liquid_saturation;
    }
    if (boost::iequals(input, "displacement"))
    {
        return Variable::displacement;
    }

    OGS_FATAL(
        "The variable name '%s' does not correspond to any known variable",
        input.c_str());

    return Variable::number_of_variables;  // to avoid the 'no return' warning
}
}  // namespace MaterialPropertyLib
