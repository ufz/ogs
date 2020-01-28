/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on January 28, 2020, 16:05 PM
 */

#include "MaterialLib/MPL/Properties/CompressibilityIdealGasLaw.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{

PropertyDataType CompressibilityIdealGasLaw::value(VariableArray const& variable_array,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t) const
{
    return 1. / std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
}

PropertyDataType CompressibilityIdealGasLaw::dValue(VariableArray const& variable_array,
                                     Variable const primary_variable,
                                     ParameterLib::SpatialPosition const& pos,
                                     double const t) const
{
    const double pressure = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);

    if (primary_variable == Variable::phase_pressure)
    {
        return -1. / (pressure * pressure);
    }

    OGS_FATAL(
        "CompressibilityIdealGasLaw::dValue is implemented for derivatives "
        "with respect to phase pressure only.");

    return 0.;
}

PropertyDataType CompressibilityIdealGasLaw::d2Value(VariableArray const& variable_array,
                                      Variable const primary_variable1,
                                      Variable const primary_variable2,
                                      ParameterLib::SpatialPosition const& pos,
                                      double const t) const
{
    const double pressure = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);

    if ((primary_variable1 == Variable::phase_pressure) &&
        (primary_variable2 == Variable::phase_pressure))
    {
        return 2. / (pressure * pressure * pressure);
    }

    OGS_FATAL(
        "CompressibilityIdealGasLaw::d2Value is implemented for derivatives "
        "with respect to phase pressure only.");

    return 0.;
}

}  // namespace MaterialPropertyLib
