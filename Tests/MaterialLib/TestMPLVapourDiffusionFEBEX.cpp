/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 5:23 PM
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Properties/VapourDiffusion/CreateVapourDiffusionFEBEX.h"
#include "TestMPL.h"

TEST(MaterialPropertyLib, VapourDiffusionFEBEX)
{
    char const xml[] =
        "<property>"
        "   <name>vapour_diffusion</name>"
        "   <type>VapourDiffusionFEBEX</type>"
        "   <tortuosity>0.8</tortuosity>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createVapourDiffusionFEBEX);
    MaterialPropertyLib::Property const& property = *property_ptr;

    double const T = 290.0;
    double const S = 0.5;
    double const phi = 0.15;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = T;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = S;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        phi;

    // The derivative of the water vapour with respect of temperature
    {
        std::array const Ts = {273.0, 293.0, 393.0, 420.0, 500.0};
        std::array const D_v_expected = {1.294719e-06, 1.470431e-06,
                                         2.494534e-06, 2.811458e-06,
                                         3.847945e-06};

        for (std::size_t i = 0; i < Ts.size(); ++i)
        {
            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = Ts[i];

            double const D_v =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(D_v_expected[i] - D_v), 1e-10)
                << "for expected water vapour diffusion " << D_v_expected[i]
                << " and for computed water vapour diffusion " << D_v;

            double const dT = 1.0e-4;
            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = Ts[i] - dT;
            double const D_v0 =
                property.template value<double>(variable_array, pos, t, dt);

            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = Ts[i] + dT;
            double const D_v1 =
                property.template value<double>(variable_array, pos, t, dt);

            double const approximated_dDv_dT = 0.5 * (D_v1 - D_v0) / dT;

            double const analytic_dDv_dT = property.template dValue<double>(
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
        std::array const D_v_expected = {
            2.886883e-06, 2.886883e-06, 2.309507e-06, 1.934212e-06,
            1.587786e-06, 1.385704e-06, 1.154753e-06, 4.330325e-07,
            0.000000e+00, 0.000000e+00};
        for (std::size_t i = 0; i < S.size(); ++i)
        {
            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = T;

            double const S_L_i = S[i];
            variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
                S_L_i;
            double const D_v =
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

            variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
                S_L_a;
            double const D_v_a =
                property.template value<double>(variable_array, pos, t, dt);
            variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
                S_L_b;
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
