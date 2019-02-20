/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VariableType.h"
#include <boost/algorithm/string/predicate.hpp>
#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
Variables convertStringToVariable(std::string const& input)
{
    if (boost::iequals(input, "phase_pressure"))
    {
        return Variables::phase_pressure;
    }
    if (boost::iequals(input, "capillary_pressure"))
    {
        return Variables::capillary_pressure;
    }
    if (boost::iequals(input, "gas_density"))
    {
        return Variables::gas_density;
    }
    if (boost::iequals(input, "liquid_density"))
    {
        return Variables::liquid_density;
    }
    if (boost::iequals(input, "temperature"))
    {
        return Variables::temperature;
    }
    if (boost::iequals(input, "liquid_saturation"))
    {
        return Variables::liquid_saturation;
    }
    if (boost::iequals(input, "u"))
    {
        return Variables::u;
    }

    OGS_FATAL(
        "The variable name '%s' does not correspond to any known variable",
        input.c_str());

    return Variables::number_of_variables;  // to avoid the 'no return' warning
}
}  // namespace MaterialPropertyLib
