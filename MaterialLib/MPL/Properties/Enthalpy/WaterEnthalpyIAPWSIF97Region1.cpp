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

#include "WaterEnthalpyIAPWSIF97Region1.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/Fluid/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{

const MaterialLib::Fluid::DimensionLessGibbsFreeEnergyRegion1
    gibbs_free_energy_;

static constexpr double ref_T_ = 1386;     ///< reference temperature in K.
static constexpr double ref_p_ = 1.653e7;  ///< reference pressure in Pa.

PropertyDataType WaterEnthalpyIAPWSIF97Region1::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = std::max(0.0, variable_array.liquid_phase_pressure);
    double const T = variable_array.temperature;

    double const tau = ref_T_ / T;
    double const pi = p / ref_p_;

    return tau *
           MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour * T *
           gibbs_free_energy_.get_dgamma_dtau(tau, pi);
}

PropertyDataType WaterEnthalpyIAPWSIF97Region1::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = std::max(0.0, variable_array.liquid_phase_pressure);
    double const T = variable_array.temperature;

    double const tau = ref_T_ / T;
    double const pi = p / ref_p_;

    if (variable == Variable::temperature)
    {
        return -tau * tau * gibbs_free_energy_.get_dgamma_dtau_dtau(tau, pi) *
               MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour;
    }

    OGS_FATAL(
        "WaterEnthalpyIAPWSIF97Region1::dValue is implemented for derivatives "
        "with respect to temperature only.");
}
}  // namespace MaterialPropertyLib
