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

#include "WaterLiquidEnthalpyIAPWSIF97Region4.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/Fluid/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"
#include "MaterialLib/MPL/Properties/WaterSaturationCurveIAPWSIF97Region4.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType WaterLiquidEnthalpyIAPWSIF97Region4::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = variable_array.liquid_phase_pressure;

    /// According to the IAPWS-IF97:
    /// http://www.iapws.org/relguide/IF97-Rev.pdf,
    /// the vapor-liquid saturation line only covers
    /// the pressure range between 611.213 Pa and 22.064 MPa.
    /// Thus, for the liquid saturation enthalpy calculated from
    /// the saturation temperature, it has the same pressure boundaries.
    if ((p < 611.213) || (p > 22.064e6))
    {
        WARN(
            "Pressure is out of the range for the water liquid saturation "
            "enthalpy.");
    }

    const MaterialLib::Fluid::DimensionLessGibbsFreeEnergyRegion1
        gibbs_free_energy_;

    static constexpr double ref_T_ = 1386;     ///< reference temperature in K.
    static constexpr double ref_p_ = 1.653e7;  ///< reference pressure in Pa.

    double const T_s = IAPWSIF97Region4::waterSaturationTemperature(p);
    double const tau = ref_T_ / T_s;
    double const pi = p / ref_p_;

    return tau *
           MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour *
           T_s * gibbs_free_energy_.get_dgamma_dtau(tau, pi);
}

PropertyDataType WaterLiquidEnthalpyIAPWSIF97Region4::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "WaterLiquidEnthalpyIAPWSIF97Region4::dValue is not implemented.");
}

}  // namespace MaterialPropertyLib
