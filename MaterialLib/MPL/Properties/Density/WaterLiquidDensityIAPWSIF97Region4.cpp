// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "WaterLiquidDensityIAPWSIF97Region4.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Properties/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"
#include "MaterialLib/MPL/Properties/WaterSaturationCurveIAPWSIF97Region4.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{

PropertyDataType WaterLiquidDensityIAPWSIF97Region4::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = variable_array.liquid_phase_pressure;

    IAPWSIF97Region4::checkPressureInRange(
        p, "the water liquid saturation density");

    const MaterialLib::Fluid::DimensionLessGibbsFreeEnergyRegion1
        gibbs_free_energy_;

    static constexpr double ref_T_ = 1386;     ///< reference temperature in K.
    static constexpr double ref_p_ = 1.653e7;  ///< reference pressure in Pa.

    double const T_s = IAPWSIF97Region4::waterSaturationTemperature(p);
    double const tau = ref_T_ / T_s;
    double const pi = p / ref_p_;

    return ref_p_ /
           (MaterialLib::PhysicalConstant::SpecificGasConstant::WaterIF97 *
            T_s * gibbs_free_energy_.get_dgamma_dpi(tau, pi));
}

PropertyDataType WaterLiquidDensityIAPWSIF97Region4::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("WaterLiquidDensityIAPWSIF97Region4::dValue is not implemented.");
}

}  // namespace MaterialPropertyLib
