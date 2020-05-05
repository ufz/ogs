/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>
#include <limits>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/GetSaturationVanGenuchten.h"
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
    double const pc_max = 20000;
    double const S_L_max = 1.0 - residual_gas_saturation;
    double const S_at_pc_max =
        std::max(residual_liquid_saturation + 1.0e-9,
                 MPL::getSaturationVanGenuchten(pc_max, p_b,
                                                residual_liquid_saturation,
                                                S_L_max, exponent));

    MPL::Property const& pressure = MPL::CapillaryPressureVanGenuchten{
        residual_liquid_saturation, residual_gas_saturation, exponent, p_b,
        S_at_pc_max};

    MPL::Property const& saturation = MPL::SaturationVanGenuchten{
        residual_liquid_saturation, residual_gas_saturation, exponent, p_b};

    MPL::VariableArray variable_array;
    fill(begin(variable_array), end(variable_array),
         std::numeric_limits<double>::quiet_NaN());
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const S_0 = -0.01;
    double const S_max = 1.01;
    int const n_steps = 10000;
    double const step_size = (S_max - S_0) / n_steps;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const S_L = S_0 + i * step_size;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L;

        double const p_cap =
            pressure.template value<double>(variable_array, pos, t, dt);
        ASSERT_LE(p_cap, pc_max);
        ASSERT_GE(p_cap, 0);

        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap;

        double const S_L_roundtrip =
            saturation.template value<double>(variable_array, pos, t, dt);

        if (S_L < residual_liquid_saturation)
        {
            EXPECT_GE(S_L_roundtrip, S_L) << "for liquid saturation " << S_L
                                          << " and capillary pressure "
                                          << p_cap;
        }

        if (residual_liquid_saturation <= S_L &&
            S_L <= 1. - residual_gas_saturation &&
            p_cap < pc_max)
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

        double const eps = 1e-8;
        double const S = std::max(S_L, S_at_pc_max);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S;

        double const dp_cap = pressure.template dValue<double>(
            variable_array, MPL::Variable::liquid_saturation, pos, t, dt);

        double const S0 =
            std::fabs(S - S_at_pc_max) < eps ? S_at_pc_max : S_L - eps;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S0;
        double const p_cap_minus =
            pressure.template value<double>(variable_array, pos, t, dt);

        double const S1 =
            std::fabs(S - S_at_pc_max) < eps ? S_at_pc_max + eps : S_L + eps;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S1;
        double const p_cap_plus =
            pressure.template value<double>(variable_array, pos, t, dt);

        double const dS = std::fabs(S - S_at_pc_max) < eps ? eps : 2.0 * eps;
        double const Dp_cap = (p_cap_plus - p_cap_minus) / dS;

        // First order derivative test

        // When p_c is close to zero, which includes the state of S >= S_max
        if (p_cap < eps)
        {
            ASSERT_LE(std::abs(dp_cap - Dp_cap), 1e-12)
                << "for liquid saturation " << S_L << ", capillary pressure "
                << p_cap << ", dp_cap/dS_L " << dp_cap << " and Dp_cap/DS_L "
                << Dp_cap;
        }
        else
        {
            if (p_cap < 0.01 * pc_max)
            {
                ASSERT_LE(std::abs((dp_cap - Dp_cap) / dp_cap), 1e-12)
                    << "for liquid saturation " << S_L
                    << ", capillary pressure " << p_cap << ", dp_cap/dS_L "
                    << dp_cap << " and Dp_cap/DS_L " << Dp_cap;
            }
            else
            {
                ASSERT_LE(std::abs((dp_cap - Dp_cap) / dp_cap), 2e-7)
                    << "for liquid saturation " << S_L
                    << ", capillary pressure " << p_cap << ", dp_cap/dS_L "
                    << dp_cap << " and Dp_cap/DS_L " << Dp_cap;
            }
        }
    }
}
