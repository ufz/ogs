/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MaterialLib/MPL/Properties/Density/CreateWaterVapourDensity.h"
#include "TestMPL.h"

static double centralDifferencesDerivative(
    double const value, double MaterialPropertyLib::VariableArray::*variable,
    double const value_increment, MaterialPropertyLib::Property const& property,
    MaterialPropertyLib::VariableArray variable_array, const auto& pos,
    double const t, double const dt)
{
    variable_array.*variable = value - value_increment;
    double const value_minus =
        property.template value<double>(variable_array, pos, t, dt);

    variable_array.*variable = value + value_increment;
    double const value_plus =
        property.template value<double>(variable_array, pos, t, dt);

    return 0.5 * (value_plus - value_minus) / value_increment;
}

TEST(MaterialPropertyLib, WaterVapourDensity)
{
    char const xml[] =
        "<property>"
        "   <name>vapour_density</name>"
        "   <type>WaterVapourDensity</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterVapourDensity);
    MaterialPropertyLib::Property const& property = *property_ptr;

    double const T = 390.0;
    double const p = -1.e+6;
    double const rho_w = 1000.0;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array.temperature = T;
    variable_array.liquid_phase_pressure = p;
    variable_array.density = rho_w;

    // The derivative of the water vapour with respect of temperature
    {
        // The first evaluation point 273.0001 is chosen to avoid lower limit
        // (T=273K) of the function in when computing numerical derivative.
        std::array const temperatures = {273.0001, 293.0, 393.0, 420.0, 500.0};
        std::array const rho_vw_expected = {4.875989e-03, 1.692871e-02,
                                            1.276865, 2.882635, 1.920386e+01};

        for (std::size_t i = 0; i < temperatures.size(); ++i)
        {
            double const T_i = temperatures[i];
            variable_array.temperature = T_i;

            double const rho_vw =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(rho_vw_expected[i] - rho_vw), 5e-6)
                << "for expected water vapour density " << rho_vw_expected[i]
                << " and for computed water vapour density " << rho_vw;

            double const approximated_drho_wv_dT = centralDifferencesDerivative(
                T_i, &MaterialPropertyLib::VariableArray::temperature, 1e-4,
                property, variable_array, pos, t, dt);

            double const analytic_drho_wv_dT = property.template dValue<double>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                t, dt);

            EXPECT_LE(std::fabs(approximated_drho_wv_dT - analytic_drho_wv_dT),
                      1e-10)
                << "for expected derivative of water vapour density with "
                   "respect to temperature "
                << approximated_drho_wv_dT
                << " and for computed derivative of water vapour density with "
                   "respect to temperature."
                << analytic_drho_wv_dT;
        }
    }

    // The derivative of the water vapour density with respect of pressure
    {
        std::array const pressures = {-1.e+7, -1.e+6, 1.e+5, 1.e+6,
                                      2.e+6,  6.e+6,  1.e+7};
        std::array const rho_vw_expected = {1.101824, 1.158320, 1.165421,
                                            1.171263, 1.177789, 1.204257,
                                            1.23132};

        for (std::size_t i = 0; i < pressures.size(); ++i)
        {
            double const p_i = pressures[i];

            variable_array.temperature = T;
            variable_array.liquid_phase_pressure = p_i;

            double const rho_vw =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(rho_vw_expected[i] - rho_vw), 5e-6)
                << "for expected water vapour density " << rho_vw_expected[i]
                << " and for computed water vapour density " << rho_vw;

            double const approximated_drho_wv_dp = centralDifferencesDerivative(
                p_i, &MaterialPropertyLib::VariableArray::liquid_phase_pressure,
                1e-3, property, variable_array, pos, t, dt);

            double const analytic_drho_wv_dp = property.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_phase_pressure, pos, t,
                dt);

            EXPECT_LE(std::fabs(approximated_drho_wv_dp - analytic_drho_wv_dp),
                      1e-13)
                << "for expected derivative of water vapour density with "
                   "respect to pressure "
                << approximated_drho_wv_dp
                << " and for computed derivative of water vapour density "
                   "with respect to pressure."
                << analytic_drho_wv_dp;
        }
    }
}
