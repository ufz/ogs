/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 19, 2021, 11:49 AM
 */

#include "WaterVapourLatentHeatWithCriticalTemperature.h"

#include <array>
#include <cmath>
#include <numeric>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
/// Critical temperature.
constexpr double T_c =
    373.92 + MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
constexpr double alpha = 1. / 8.;
constexpr double beta = 1. / 3.;
constexpr double Delta = 0.79 - beta;
// Coefficients a_1, a_2, a_4, b_1, ... b_5
constexpr std::array c = {1989.41582,    11178.45586,  26923.68994,
                          -28989.28947,  -19797.03646, 28403.32283,
                          -30382.306422, 15210.380};

PropertyDataType WaterVapourLatentHeatWithCriticalTemperature::value(
    const VariableArray& variable_array,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*t*/,
    const double /*dt*/) const
{
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    if (T >= T_c)
    {
        return 0.0;
    }

    double const tau = (T_c - T) / T_c;
    double const tau_p2 = tau * tau;
    double const tau_p3 = tau_p2 * tau;
    double const tau_p4 = tau_p3 * tau;
    double const tau_p5 = tau_p4 * tau;

    std::array v = {std::pow(tau, beta),
                    std::pow(tau, beta + Delta),
                    std::pow(tau, 1 - alpha + beta),
                    tau,
                    tau_p2,
                    tau_p3,
                    tau_p4,
                    tau_p5};

    // The formula gives the value in kJ/kg, and the return value is in the
    // units of J/kg.

    return 1000.0 *
#if __GNUC__ < 9 || (__GNUC__ == 9 && (__GNUC_MINOR__ < 3))
           std::inner_product
#else
           std::transform_reduce
#endif
           (begin(c), end(c), begin(v), 0.);
}

PropertyDataType WaterVapourLatentHeatWithCriticalTemperature::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (primary_variable != Variable::temperature)
    {
        OGS_FATAL(
            "WaterVapourLatentHeatWithCriticalTemperature::dValue is "
            "implemented "
            "for "
            "the derivative with respect to temperature only.");
    }
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    if (T >= T_c)
    {
        return 0.0;
    }

    constexpr std::array dc = {c[0] * beta,
                               c[1] * (beta + Delta),
                               c[2] * (1 - alpha + beta),
                               c[3],
                               c[4] * 2,
                               c[5] * 3,
                               c[6] * 4,
                               c[7] * 5};

    double const tau = (T_c - T) / T_c;
    double const tau_p2 = tau * tau;
    double const tau_p3 = tau_p2 * tau;
    double const tau_p4 = tau_p3 * tau;
    std::array v = {std::pow(tau, beta - 1),
                    std::pow(tau, beta + Delta - 1),
                    std::pow(tau, -alpha + beta),
                    1.,
                    tau,
                    tau_p2,
                    tau_p3,
                    tau_p4};

    // The formula gives the value in kJ/kg/K, and the value is return in
    // the unit of J/kg/K.
    return -1000.0 *
#if __GNUC__ < 9 || (__GNUC__ == 9 && (__GNUC_MINOR__ < 3))
           std::inner_product
#else
           std::transform_reduce
#endif
           (begin(dc), end(dc), begin(v), 0.) /
           T_c;
}

}  // namespace MaterialPropertyLib
