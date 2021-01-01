/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationVanGenuchten.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, CapillaryPressureVanGenuchten)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.79;
    double const p_b = 10000;
    double const maximum_capillary_pressure = 20000;

    MPL::Property const& pressure =
        MPL::CapillaryPressureVanGenuchten{"capillary_pressure",
                                           residual_liquid_saturation,
                                           residual_gas_saturation,
                                           exponent,
                                           p_b,
                                           maximum_capillary_pressure};

    MPL::Property const& saturation =
        MPL::SaturationVanGenuchten{"saturation", residual_liquid_saturation,
                                    residual_gas_saturation, exponent, p_b};

    MPL::VariableArray variable_array;
    fill(begin(variable_array), end(variable_array),
         std::numeric_limits<double>::quiet_NaN());
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const S_0 = -0.01;
    double const S_max = 1.01;
    int const n_steps = 10000;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const S_L = S_0 + i * (S_max - S_0) / n_steps;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L;

        double const p_cap =
            pressure.template value<double>(variable_array, pos, t, dt);
        ASSERT_LE(p_cap, maximum_capillary_pressure);
        ASSERT_GE(p_cap, 0);

        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap;

        double const S_L_roundtrip =
            saturation.template value<double>(variable_array, pos, t, dt);

        if (S_L < residual_liquid_saturation)
        {
            EXPECT_GE(S_L_roundtrip, S_L)
                << "for liquid saturation " << S_L << " and capillary pressure "
                << p_cap;
        }

        if (residual_liquid_saturation <= S_L &&
            S_L <= 1. - residual_gas_saturation &&
            p_cap < maximum_capillary_pressure)
        {
            ASSERT_LE(std::abs(S_L - S_L_roundtrip), 1e-15)
                << "for liquid saturation " << S_L << " and capillary pressure "
                << p_cap;
        }

        if (1. - residual_gas_saturation < S_L)
        {
            ASSERT_EQ(S_L_roundtrip, 1. - residual_gas_saturation)
                << "for liquid saturation " << S_L << " and capillary pressure "
                << p_cap;
        }

        double const dp_cap = pressure.template dValue<double>(
            variable_array, MPL::Variable::liquid_saturation, pos, t, dt);

        double const eps = 1e-6;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L - eps;
        double const p_cap_minus =
            pressure.template value<double>(variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L + eps;
        double const p_cap_plus =
            pressure.template value<double>(variable_array, pos, t, dt);

        double const Dp_cap = (p_cap_plus - p_cap_minus) / 2 / eps;

        //
        // First order derivative tests
        //

        // Testing relative error up to the S_L_max with different tolerances as
        // dp_cap/dS_L approaches -infinity.
        if (S_L < 1 - 2 * residual_gas_saturation)
        {
            ASSERT_LE(std::abs(dp_cap - Dp_cap) / p_cap, 1e-9)
                << "for liquid saturation " << S_L << ", capillary pressure "
                << p_cap << ", dp_cap/dS_L " << dp_cap << " and Dp_cap/DS_L "
                << Dp_cap;
        }
        else if (S_L < 1 - residual_gas_saturation * 1.35)
        {
            ASSERT_LE(std::abs(dp_cap - Dp_cap) / p_cap, 1e-8)
                << "for liquid saturation " << S_L << ", capillary pressure "
                << p_cap << ", dp_cap/dS_L " << dp_cap << " and Dp_cap/DS_L "
                << Dp_cap;
        }
        // Skip the range up to maximum liquid saturation, continue after that
        // with absolute error.
        else if (S_L > 1 - residual_gas_saturation)
        {
            ASSERT_EQ(dp_cap, Dp_cap)
                << "for liquid saturation " << S_L << ", capillary pressure "
                << p_cap << ", dp_cap/dS_L " << dp_cap << " and Dp_cap/DS_L "
                << Dp_cap;
        }
    }
}
