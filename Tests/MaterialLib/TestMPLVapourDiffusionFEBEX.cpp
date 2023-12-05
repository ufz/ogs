/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
        "   <base_diffusion_coefficient>2.16e-5</base_diffusion_coefficient>"
        "   <exponent>1.8</exponent>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createVapourDiffusionFEBEX);
    MaterialPropertyLib::Property const& property = *property_ptr;

    double const T = 290.0;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array.temperature = T;

    // The derivative of the water vapour with respect of temperature
    {
        std::array const Ts = {273.0, 293.0, 393.0, 420.0, 500.0};
        std::array const D_v_expected = {1.72629e-05 / 0.8, 1.96057e-05 / 0.8,
                                         3.32605e-05 / 0.8, 3.74861e-05 / 0.8,
                                         5.13059e-05 / 0.8, 1.92459e-05 / 0.8};

        for (std::size_t i = 0; i < Ts.size(); ++i)
        {
            variable_array.temperature = Ts[i];

            double const D_v =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(D_v_expected[i] - D_v), 1e-10)
                << "for expected water vapour diffusion " << D_v_expected[i]
                << " and for computed water vapour diffusion " << D_v;

            double const dT = 1.0e-4;
            variable_array.temperature = Ts[i] - dT;
            double const D_v0 =
                property.template value<double>(variable_array, pos, t, dt);

            variable_array.temperature = Ts[i] + dT;
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
    // The derivative of the water vapour with respect of saturation, which is
    // zero.
    {
        std::array const S = {-1.0, 0.0, 0.2,  0.33, 0.45,
                              0.52, 0.6, 0.85, 1.0,  1.1};
        double const D_v_expected = 1.92459e-05 / 0.8;
        for (std::size_t i = 0; i < S.size(); ++i)
        {
            variable_array.temperature = T;

            double const D_v =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(D_v_expected - D_v), 1e-10)
                << "for expected water vapour diffusion " << D_v_expected
                << " and for computed water vapour diffusion " << D_v;

            double const analytic_dDv_dS = property.template dValue<double>(
                variable_array, MPL::Variable::liquid_saturation, pos, t, dt);

            double const expected_dDv_dS = 0.0;

            ASSERT_LE(std::fabs(expected_dDv_dS - analytic_dDv_dS), 1e-10)
                << "for expected derivative of water vapour diffusion with "
                   "respect to saturation "
                << expected_dDv_dS
                << " and for computed derivative of water vapour diffusion "
                   "with respect to saturation."
                << analytic_dDv_dS;
        }
    }
}
