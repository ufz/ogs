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

#include "WaterLiquidEnthalpyIAPWSIF97Region4.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/WaterSaturationCurveIAPWSIF97Region4.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType WaterLiquidEnthalpyIAPWSIF97Region4::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = variable_array.phase_pressure;
    if (p < 0)
    {
        OGS_FATAL(
            "WaterLiquidEnthalpyIAPWSIF97Region4 can not be calculated from "
            "a negative phase pressure value.");
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
