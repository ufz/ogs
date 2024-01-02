/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
    double const p = variable_array.liquid_phase_pressure;

    /// According to the IAPWS-IF97:
    /// http://www.iapws.org/relguide/IF97-Rev.pdf,
    /// the vapor-liquid saturation line in region4 only covers
    /// the pressure range between 611.213 Pa and 22.064 MPa.
    if ((p < 611.213) || (p > 22.064e6))
    {
        WARN(
            "Pressure is out of the range for the water saturation temperature "
            "in region4.");
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
