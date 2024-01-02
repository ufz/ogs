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

#include "WaterTemperatureIAPWSIF97Region1.h"

#include <cmath>

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
static constexpr std::array n_T = {
    -0.23872489924521e3,  0.40421188637945e3,   0.11349746881718e3,
    -0.58457616048039e1,  -0.15285482413140e-3, -0.10866707695377e-5,
    -0.13391744872602e2,  0.43211039183559e2,   -0.54010067170506e2,
    0.30535892203916e2,   -0.65964749423638e1,  0.93965400878363e-2,
    0.11573647505340e-6,  -0.25858641282073e-4, -0.40644363084799e-8,
    0.66456186191635e-7,  0.80670734103027e-10, -0.93477771213947e-12,
    0.58265442020601e-14, -0.15020185953503e-16};

static constexpr std::array I_T = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                                   1, 1, 1, 2, 2, 3, 3, 4, 5, 6};

static constexpr std::array J_T = {0, 1,  2,  6,  22, 32, 0,  1,  2,  3,
                                   4, 10, 32, 10, 32, 10, 32, 32, 32, 32};

double computeTemperature(double const pi, double const eta)
{
    double val = 0.;
    for (int i = 0; i < 20; i++)
    {
        val += n_T[i] * std::pow(pi, I_T[i]) * std::pow(eta + 1, J_T[i]);
    }

    return val;
}

PropertyDataType WaterTemperatureIAPWSIF97Region1::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p = std::max(0.0, variable_array.liquid_phase_pressure);
    double const h = variable_array.enthalpy;

    static constexpr double ref_h_ = 2500e3;  ///< reference enthalpy in J/kg.
    static constexpr double ref_p_ = 1e6;     ///< reference pressure in Pa.

    double const eta = h / ref_h_;
    double const pi = p / ref_p_;

    return computeTemperature(pi, eta);
}

PropertyDataType WaterTemperatureIAPWSIF97Region1::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("WaterTemperatureIAPWSIF97Region1::dValue is not implemented.");
}

}  // namespace MaterialPropertyLib
