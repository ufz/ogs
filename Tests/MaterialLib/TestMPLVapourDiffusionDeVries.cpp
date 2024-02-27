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

#include "MaterialLib/MPL/Properties/VapourDiffusion/CreateVapourDiffusionDeVries.h"
#include "TestMPL.h"

TEST(MaterialPropertyLib, VapourDiffusionDeVries)
{
    char const xml[] =
        "<property>"
        "   <name>vapour_diffusion</name>"
        "   <type>VapourDiffusionDeVries</type>"
        "   <base_diffusion_coefficient>2.16e-5</base_diffusion_coefficient>"
        "   <exponent>1.8</exponent>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createVapourDiffusionDeVries);
    MaterialPropertyLib::Property const& property = *property_ptr;

    double const T = 290.0;
    double const pG = 1.0e5;
    double const S = 0.3;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array.temperature = T;
    variable_array.gas_phase_pressure = pG;
    variable_array.liquid_saturation = S;

    // The derivative of the water vapour with respect of saturation, which is
    // zero.
    double const S_0 = -0.1;
    double const S_max = 1.1;
    int n_steps = 10;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const S_i = S_0 + i * (S_max - S_0) / n_steps;
        variable_array.liquid_saturation = S_i;

        double const dD_v = property.template dValue<double>(
            variable_array, MPL::Variable::liquid_saturation, pos, t, dt);

        double const eps = 1e-8;
        variable_array.liquid_saturation = S_i - eps;
        double const D_v_minus =
            property.template value<double>(variable_array, pos, t, dt);
        variable_array.liquid_saturation = S_i + eps;
        double const D_v_plus =
            property.template value<double>(variable_array, pos, t, dt);

        double const DD_v = (D_v_plus - D_v_minus) / 2 / eps;

        ASSERT_LE(std::fabs(dD_v - DD_v), 0.0)
            << "for expected derivative of water vapour diffusion with "
               "respect to saturation "
            << DD_v
            << " and for computed derivative of water vapour diffusion "
               "with respect to saturation."
            << dD_v;
    }
    // The derivative of the water vapour with respect of temperature
    double const T_0 = 273.15;
    double const T_max = 500;
    n_steps = 1000;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const T_i = T_0 + i * (T_max - S_0) / n_steps;
        variable_array.temperature = T_i;

        double const dD_v = property.template dValue<double>(
            variable_array, MPL::Variable::temperature, pos, t, dt);

        double const eps = 1e-8;
        variable_array.temperature = T_i - eps;
        double const D_v_minus =
            property.template value<double>(variable_array, pos, t, dt);
        variable_array.temperature = T_i + eps;
        double const D_v_plus =
            property.template value<double>(variable_array, pos, t, dt);

        double const DD_v = (D_v_plus - D_v_minus) / 2 / eps;

        ASSERT_LE(std::fabs(dD_v - DD_v), 1e-9)
            << "for expected derivative of water vapour diffusion with "
               "respect to temperature "
            << DD_v
            << " and for computed derivative of water vapour diffusion "
               "with respect to temperature."
            << dD_v;
    }
    // The derivative of the water vapour with respect of gas pressure
    double const pG_0 = 1.0;
    double const pG_max = 1.0e7;
    n_steps = 10000;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const pG_i = pG_0 + i * (pG_max - pG_0) / n_steps;
        variable_array.gas_phase_pressure = pG_i;

        double const dD_v = property.template dValue<double>(
            variable_array, MPL::Variable::gas_phase_pressure, pos, t, dt);

        double const eps = 1e-8;
        variable_array.gas_phase_pressure = pG_i - eps;
        double const D_v_minus =
            property.template value<double>(variable_array, pos, t, dt);
        variable_array.gas_phase_pressure = pG_i + eps;
        double const D_v_plus =
            property.template value<double>(variable_array, pos, t, dt);

        double const DD_v = (D_v_plus - D_v_minus) / 2 / eps;

        ASSERT_LE(std::fabs(dD_v - DD_v), 5e-9)
            << "for expected derivative of water vapour diffusion with "
               "respect to gas pressure "
            << DD_v
            << " and for computed derivative of water vapour diffusion "
               "with respect to gas pressure."
            << dD_v;
    }
}
