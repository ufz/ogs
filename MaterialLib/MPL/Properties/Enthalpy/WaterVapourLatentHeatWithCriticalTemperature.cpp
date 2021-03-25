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

#include <algorithm>
#include <array>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
static std::array<double, 3> constexpr a = {1989.41582, 11178.45586,
                                            26923.68994};
static std::array<double, 3> constexpr a_exp = {0.3333333333333333, 0.79,
                                                1.2083333333333333};

static std::array<double, 5> constexpr b = {
    -28989.28947, -19797.03646, 28403.32283, -30382.306422, 15210.380};

PropertyDataType WaterVapourLatentHeatWithCriticalTemperature::value(
    const VariableArray& variable_array,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*t*/,
    const double /*dt*/) const
{
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    if (T >= T_c_)
    {
        return 0.0;
    }

    double const tau = (T_c_ - T) / T_c_;
    double L_w = 0.0;
    for (std::size_t i = 0; i < a.size(); i++)
    {
        L_w += a[i] * std::pow(tau, a_exp[i]);
    }

    double tau_to_power_i = tau;

    std::for_each(b.begin(), b.end(), [&](double const b_i) {
        L_w += b_i * tau_to_power_i;  // b_i * tau^i
        tau_to_power_i *= tau;
    });

    // The formula gives the value in kJ/kg, and the value is return in
    // the unit of J/kg.
    return 1000.0 * L_w;
}

PropertyDataType WaterVapourLatentHeatWithCriticalTemperature::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (primary_variable == Variable::temperature)
    {
        const double T = std::get<double>(
            variable_array[static_cast<int>(Variable::temperature)]);

        if (T >= T_c_)
        {
            return 0.0;
        }

        double const tau = (T_c_ - T) / T_c_;
        double dL_w_dtau = 0.0;
        for (std::size_t i = 0; i < a.size(); i++)
        {
            dL_w_dtau += a[i] * a_exp[i] * std::pow(tau, a_exp[i] - 1.0);
        }

        double tau_to_power_i = 1.0;
        int exponent_index = 0;
        std::for_each(b.begin(), b.end(), [&](double const b_i) {
            dL_w_dtau +=
                (exponent_index + 1) * b_i * tau_to_power_i;  // b_i * dau^i/dT
            tau_to_power_i *= tau;
            exponent_index++;
        });

        // The formula gives the value in kJ/kg/K, and the value is return in
        // the unit of J/kg/K.
        return -1000.0 * dL_w_dtau / T_c_;
    }

    OGS_FATAL(
        "WaterVapourLatentHeatWithCriticalTemperature::dValue is implemented "
        "for "
        "the derivative with respect to temperature only.");
}

}  // namespace MaterialPropertyLib
