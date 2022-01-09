/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationBrooksCorey.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, SaturationBrooksCorey)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 3.0;
    double const p_b = 5000;
    MPL::Property const& pressure_saturation =
        MPL::SaturationBrooksCorey{"saturation", residual_liquid_saturation,
                                   residual_gas_saturation, exponent, p_b};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const p_0 = -1e6;
    double const p_max = 4000;
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

        double const eps = 1e-2;
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
