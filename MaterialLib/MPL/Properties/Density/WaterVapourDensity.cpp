/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#include "WaterVapourDensity.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
/// \f$\rho_{vS}\f$
static double saturatedVaporDensity(double const T)
{
    return 1.0e-3 * std::exp(19.819 - 4975.9 / T);
}

/// \f$\frac{\partial \rho_{vS}}{\partial T}\f$
static double dsaturatedVaporDensitydT(double const T)
{
    return 4.9759 * std::exp(19.819 - 4975.9 / T) / (T * T);
}

static double humidity(double const T, double const p,
                       double const water_density)
{
    return std::exp(
        p / (MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour *
             T * water_density));
}

PropertyDataType WaterVapourDensity::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    double const T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    double const water_density =
        std::get<double>(variable_array[static_cast<int>(Variable::density)]);

    return humidity(T, p, water_density) * saturatedVaporDensity(T);
}

PropertyDataType WaterVapourDensity::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    double const T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    double const water_density =
        std::get<double>(variable_array[static_cast<int>(Variable::density)]);

    if (primary_variable == Variable::temperature)
    {
        double const h = humidity(T, p, water_density);
        double const rho_v = h * saturatedVaporDensity(T);
        double const drho_vS_dT = dsaturatedVaporDensitydT(T);

        return h * drho_vS_dT - rho_v * p /
                                    (water_density * T * T *
                                     MaterialLib::PhysicalConstant::
                                         SpecificGasConstant::WaterVapour);
    }

    if (primary_variable == Variable::phase_pressure)
    {
        double const factor =
            1.0 /
            (MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour *
             T * water_density);
        double const dhumiditydp = factor * std::exp(factor * p);

        return dhumiditydp * saturatedVaporDensity(T);
    }

    OGS_FATAL(
        "WaterVapourDensity::dValue is implemented for derivatives with "
        "respect to temperature or phase_pressure only.");
}

}  // namespace MaterialPropertyLib
