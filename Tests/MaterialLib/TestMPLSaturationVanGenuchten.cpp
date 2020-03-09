/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <random>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/SaturationVanGenuchten.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, SaturationVanGenuchten)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.79;
    double const entry_pressure = 5000;
    double const max_capillary_pressure = std::numeric_limits<double>::max();

    MPL::Property const& pressure_saturation = MPL::SaturationVanGenuchten{
        residual_liquid_saturation, residual_gas_saturation, exponent,
        entry_pressure, max_capillary_pressure};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const p_0 = -1e6;
    double const p_max = 10000;
    int const n_steps = 10000;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const p_L = p_0 + i * (p_max - p_0) / n_steps;
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            -p_L;

        double const S = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);
        double const dS = pressure_saturation.template dValue<double>(
            variable_array, MPL::Variable::capillary_pressure, pos, t, dt);
        double const dS2 = pressure_saturation.template d2Value<double>(
            variable_array, MPL::Variable::capillary_pressure,
            MPL::Variable::capillary_pressure, pos, t, dt);

        double const eps = 1e-1;
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            -p_L - eps;
        double const S_minus = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            -p_L + eps;
        double const S_plus = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);

        double const DS = (S_plus - S_minus) / 2 / eps;
        double const DS2 = (S_plus - 2 * S + S_minus) / (eps * eps);

        ASSERT_LE(std::abs(dS - DS), 1e-9)
            << "for capillary pressure " << -p_L << " and saturation " << S;
        ASSERT_LE(std::abs(dS2 - DS2), 1e-9)
            << "for capillary pressure " << -p_L << " and saturation " << S;
    }
}

TEST(MaterialPropertyLib, CapillaryPressureVanGenuchten)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.79;
    double const entry_pressure = 5000;
    double const max_capillary_pressure = std::numeric_limits<double>::max();

    MPL::Property const& pressure_saturation = MPL::SaturationVanGenuchten{
        residual_liquid_saturation, residual_gas_saturation, exponent,
        entry_pressure, max_capillary_pressure};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::random_device rd;
    std::mt19937 mt(rd());
    const double offset = std::sqrt(std::numeric_limits<double>::epsilon());
    std::uniform_real_distribution<double> distributor(
        residual_liquid_saturation + offset,
        1.0 - residual_gas_saturation - offset);

    const int n = 20;
    for (int i = 0; i <= n; ++i)
    {
        double const S = distributor(mt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S;

        const double computed_capillary_pressure =
            pressure_saturation.template inverse_value<double>(variable_array,
                                                               pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            computed_capillary_pressure;

        const double re_computedS = pressure_saturation.template value<double>(
            variable_array, pos, t, dt);

        ASSERT_LE(std::abs(S - re_computedS), 1e-9)
            << "for saturation " << S
            << " and re-computed saturation via capillary pressure"
            << re_computedS;
    }
}
