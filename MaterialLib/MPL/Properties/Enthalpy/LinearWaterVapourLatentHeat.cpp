/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 16, 2021, 10:03 AM
 */

#include "LinearWaterVapourLatentHeat.h"

#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType LinearWaterVapourLatentHeat::value(
    const VariableArray& variable_array,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*t*/,
    const double /*dt*/) const
{
    const double T = variable_array.temperature;

    return 2.501e+6 -
           2369.2 * (T - MaterialLib::PhysicalConstant::CelsiusZeroInKelvin);
}

PropertyDataType LinearWaterVapourLatentHeat::dValue(
    VariableArray const& /*variable_array*/, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable == Variable::temperature)
    {
        return -2369.2;
    }

    OGS_FATAL(
        "LinearWaterVapourLatentHeat::dValue is implemented for the derivative "
        "with respect to temperature only.");
}

}  // namespace MaterialPropertyLib
