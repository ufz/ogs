// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "WaterVapourDensityIAPWSIF97Region4.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Properties/GibbsFreeEnergy/DimensionlessGibbsFreeEnergyRegion2.h"
#include "MaterialLib/MPL/Properties/WaterSaturationCurveIAPWSIF97Region4.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{

PropertyDataType WaterVapourDensityIAPWSIF97Region4::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = variable_array.liquid_phase_pressure;

    IAPWSIF97Region4::checkPressureInRange(
        p, "the water vapour saturation density");

    double const T_s = IAPWSIF97Region4::waterSaturationTemperature(p);

    static constexpr double ref_T_ = 540;   ///< reference temperature in K.
    static constexpr double ref_p_ = 1.e6;  ///< reference pressure in Pa.
    double const tau = ref_T_ / T_s;
    double const pi = p / ref_p_;

    double dgamma_dtau =
        MaterialLib::Fluid::DimensionlessGibbsFreeEnergyRegion2::getdGammadPi(
            tau, pi);

    return p /
           (pi * MaterialLib::PhysicalConstant::SpecificGasConstant::WaterIF97 *
            T_s * dgamma_dtau);
}

PropertyDataType WaterVapourDensityIAPWSIF97Region4::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("WaterVapourDensityIAPWSIF97Region4::dValue is not implemented.");
}

}  // namespace MaterialPropertyLib
