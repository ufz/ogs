/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#include "WaterDensityIAPWSIF97Region1.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{

PropertyDataType WaterDensityIAPWSIF97Region1::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = std::max(0.0, variable_array.liquid_phase_pressure);
    double const T = variable_array.temperature;
    const double tau = ref_T_ / T;
    const double pi = p / ref_p_;

    return ref_p_ /
           (MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour *
            T * gibbs_free_energy_.get_dgamma_dpi(tau, pi));
}

PropertyDataType WaterDensityIAPWSIF97Region1::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double T = variable_array.temperature;
    double const p = std::max(0.0, variable_array.liquid_phase_pressure);

    const double tau = ref_T_ / T;
    const double pi = p / ref_p_;

    const double dgamma_dpi = gibbs_free_energy_.get_dgamma_dpi(tau, pi);

    switch (variable)
    {
        case Variable::temperature:
            return -(ref_p_ -
                     tau * ref_p_ *
                         gibbs_free_energy_.get_dgamma_dtau_dpi(tau, pi) /
                         dgamma_dpi) /
                   (MaterialLib::PhysicalConstant::SpecificGasConstant::
                        WaterVapour *
                    T * T * dgamma_dpi);
        case Variable::liquid_phase_pressure:
            return -gibbs_free_energy_.get_dgamma_dpi_dpi(tau, pi) /
                   (MaterialLib::PhysicalConstant::SpecificGasConstant::
                        WaterVapour *
                    T * dgamma_dpi * dgamma_dpi);
        default:
            OGS_FATAL(
                "WaterDensityIAPWSIF97Region1::dValue is implemented for "
                "derivatives with "
                "respect to temperature or liquid_phase_pressure only.");
    }
}

}  // namespace MaterialPropertyLib
