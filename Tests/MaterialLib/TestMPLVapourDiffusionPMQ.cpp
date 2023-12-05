/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 7, 2021, 9:20 AM
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Properties/VapourDiffusion/CreateVapourDiffusionPMQ.h"
#include "TestMPL.h"

TEST(MaterialPropertyLib, VapourDiffusionPMQ)
{
    const char xml[] =
        "<property>"
        "   <name>vapour_diffusion</name>"
        "   <type>VapourDiffusionPMQ</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createVapourDiffusionPMQ);
    MaterialPropertyLib::Property const& property = *property_ptr;

    const double T = 290.0;
    const double S = 0.5;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array.temperature = T;
    variable_array.liquid_saturation = S;

    // The derivative of the water vapour with respect of temperature
    {
        std::array const Ts = {273.0, 293.0, 393.0, 420.0, 500.0};
        std::array const D_v_expected = {
            7.1209560000000006e-06 / 0.66, 8.087368e-06 / 0.66,
            1.3719933333333336e-05 / 0.66, 1.5463013333333333e-05 / 0.66,
            2.1163693333333336e-05 / 0.66};

        for (std::size_t i = 0; i < Ts.size(); ++i)
        {
            variable_array.temperature = Ts[i];

            const double D_v =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(D_v_expected[i] - D_v), 1e-10)
                << "for expected water vapour diffusion " << D_v_expected[i]
                << " and for computed water vapour diffusion " << D_v;

            const double dT = 1.0e-4;
            variable_array.temperature = Ts[i] - dT;
            const double D_v0 =
                property.template value<double>(variable_array, pos, t, dt);

            variable_array.temperature = Ts[i] + dT;
            const double D_v1 =
                property.template value<double>(variable_array, pos, t, dt);

            const double approximated_dDv_dT = 0.5 * (D_v1 - D_v0) / dT;

            const double analytic_dDv_dT = property.template dValue<double>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                t, dt);

            ASSERT_LE(std::fabs(approximated_dDv_dT - analytic_dDv_dT), 1e-7)
                << "for expected derivative of water vapour diffusion with "
                   "respect to temperature "
                << approximated_dDv_dT
                << " and for computed derivative of water vapour diffusion "
                   "with respect to temperature."
                << analytic_dDv_dT;
        }
    }
    // The derivative of the water vapour with respect of saturation
    {
        std::array const S = {-1.0, 0.0, 0.2,  0.33, 0.45,
                              0.52, 0.6, 0.85, 1.0,  1.1};
        std::array const D_v_expected = {1.58779e-05 / 0.66,
                                         1.58779e-05 / 0.66,
                                         1.27023e-05 / 0.66,
                                         1.06382e-05 / 0.66,
                                         8.73282e-06 / 0.66,
                                         7.62137e-06 / 0.66,
                                         6.35114e-06 / 0.66,
                                         2.38168e-06 / 0.66,
                                         0.0,
                                         0.0};
        for (std::size_t i = 0; i < S.size(); ++i)
        {
            variable_array.temperature = T;

            double const S_L_i = S[i];
            variable_array.liquid_saturation = S_L_i;
            const double D_v =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(D_v_expected[i] - D_v), 1e-10)
                << "for expected water vapour diffusion " << D_v_expected[i]
                << " and for computed water vapour diffusion " << D_v;

            double const analytic_dDv_dS = property.template dValue<double>(
                variable_array, MPL::Variable::liquid_saturation, pos, t, dt);
            double const dS_L = 1e-8;
            double S_L_a = S_L_i - dS_L;
            double S_L_b = S_L_i + dS_L;
            double factor = 0.5;

            if (S_L_i <= 0.0)
            {
                S_L_a = S_L_i;
                S_L_b = dS_L;
                factor = 1.0;
            }
            else if (S_L_i >= 1.0)
            {
                S_L_a = 1.0 - dS_L;
                S_L_b = S_L_i;
                factor = 1.0;
            }

            variable_array.liquid_saturation = S_L_a;
            double const D_v_a =
                property.template value<double>(variable_array, pos, t, dt);
            variable_array.liquid_saturation = S_L_b;
            double const D_v_b =
                property.template value<double>(variable_array, pos, t, dt);

            double const approximated_dDv_dS = factor * (D_v_b - D_v_a) / dS_L;

            ASSERT_LE(std::fabs(approximated_dDv_dS - analytic_dDv_dS) /
                          analytic_dDv_dS,
                      1e-10)
                << "for expected derivative of water vapour diffusion with "
                   "respect to saturation "
                << approximated_dDv_dS
                << " and for computed derivative of water vapour diffusion "
                   "with respect to saturation."
                << analytic_dDv_dS;
        }
    }
}
