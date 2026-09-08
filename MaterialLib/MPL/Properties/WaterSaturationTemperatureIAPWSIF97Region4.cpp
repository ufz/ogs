// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

    IAPWSIF97Region4::checkPressureInRange(
        p, "the water saturation temperature in region4");

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
