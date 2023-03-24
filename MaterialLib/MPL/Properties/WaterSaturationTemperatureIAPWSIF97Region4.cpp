/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on Feb 8, 2023, 3:05 PM
 */

#include "WaterSaturationTemperatureIAPWSIF97Region4.h"

#include "BaseLib/Error.h"
#include "WaterSaturationCurveIAPWSIF97Region4.h"

namespace MaterialPropertyLib
{
PropertyDataType WaterSaturationTemperatureIAPWSIF97Region4::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = variable_array.phase_pressure;
    if (p < 0)
    {
        OGS_FATAL("The phase pressure can not be a negative value.");
    }

    return IAPWSIF97Region4::waterSaturationTemperature(p);
}

PropertyDataType WaterSaturationTemperatureIAPWSIF97Region4::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "WaterSaturationTemperatureIAPWSIF97Region4::dValue is not "
        "implemented.");
}

}  // namespace MaterialPropertyLib
