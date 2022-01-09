/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = T;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)] = p;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::density)] =
        rho_w;

    // The derivative of the water vapour with respect of temperature
    {
        std::array const temperatures = {273.0, 293.0, 393.0, 420.0, 500.0};
        std::array const rho_vw_expected = {4.875989e-03, 1.692871e-02,
                                            1.276865, 2.882635, 1.920386e+01};

        for (std::size_t i = 0; i < temperatures.size(); ++i)
        {
            double const T_i = temperatures[i];
            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = T_i;

            double const rho_vw =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(rho_vw_expected[i] - rho_vw), 5e-6)
                << "for expected water vapour density " << rho_vw_expected[i]
                << " and for computed water vapour density " << rho_vw;

            double const dT = 1.0e-4;
            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = T_i - dT;
            double const rho_vw0 =
                property.template value<double>(variable_array, pos, t, dt);

            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = T_i + dT;
            double const rho_vw1 =
                property.template value<double>(variable_array, pos, t, dt);

            double const approximated_drho_wv_dT =
                0.5 * (rho_vw1 - rho_vw0) / dT;

            double const analytic_drho_wv_dT = property.template dValue<double>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                t, dt);

            ASSERT_LE(std::fabs(approximated_drho_wv_dT - analytic_drho_wv_dT),
                      1e-6)
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

            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::temperature)] = T;
            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_i;

            double const rho_vw =
                property.template value<double>(variable_array, pos, t, dt);

            ASSERT_LE(std::fabs(rho_vw_expected[i] - rho_vw), 5e-6)
                << "for expected water vapour density " << rho_vw_expected[i]
                << " and for computed water vapour density " << rho_vw;

            double const dp = 1.0e-3;
            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_i - dp;

            double const rho_vw0 =
                property.template value<double>(variable_array, pos, t, dt);

            variable_array[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_i + dp;
            double const rho_vw1 =
                property.template value<double>(variable_array, pos, t, dt);

            double const approximated_drho_wv_dp =
                0.5 * (rho_vw1 - rho_vw0) / dp;

            double const analytic_drho_wv_dp = property.template dValue<double>(
                variable_array, MaterialPropertyLib::Variable::phase_pressure,
                pos, t, dt);

            ASSERT_LE(std::fabs(approximated_drho_wv_dp - analytic_drho_wv_dp),
                      1e-7)
                << "for expected derivative of water vapour density with "
                   "respect to pressure "
                << approximated_drho_wv_dp
                << " and for computed derivative of water vapour density "
                   "with respect to pressure."
                << analytic_drho_wv_dp;
        }
    }
}
