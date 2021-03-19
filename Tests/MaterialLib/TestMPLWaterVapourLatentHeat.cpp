/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 16, 2021, 10:44 AM
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Properties/Enthalpy/CreateLinearWaterVapourLatentHeat.h"
#include "TestMPL.h"

TEST(MaterialPropertyLib, LinearWaterVapourLatentHeat)
{
    char const xml[] =
        "<property>"
        "   <name>latent_heat</name>"
        "   <type>LinearWaterVapourLatentHeat</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createLinearWaterVapourLatentHeat);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::array const temperatures = {273.0, 293.0, 393.0, 420.0, 500.0};
    std::array const L_wv_expected = {2.501355e+06, 2.453971e+06, 2.217051e+06,
                                      2.153083e+06, 1.963547e+06};

    for (std::size_t i = 0; i < temperatures.size(); ++i)
    {
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = temperatures[i];

        double const L_wv =
            property.template value<double>(variable_array, pos, t, dt);

        ASSERT_LE(std::fabs(L_wv_expected[i] - L_wv) / L_wv_expected[i], 1e-6)
            << "for expected water vapour latent heat " << L_wv_expected[i]
            << " and for computed water vapour latent heat " << L_wv;

        double const dT = 1.0e-4;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = temperatures[i] - dT;
        double const L_wv0 =
            property.template value<double>(variable_array, pos, t, dt);

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = temperatures[i] + dT;
        double const L_wv1 =
            property.template value<double>(variable_array, pos, t, dt);

        double const approximated_dLw_dT = 0.5 * (L_wv1 - L_wv0) / dT;

        double const analytic_dLw_dT = property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::temperature, pos, t,
            dt);

        ASSERT_LE(std::fabs(approximated_dLw_dT - analytic_dLw_dT), 1e-6)
            << "for expected derivative of water vapour latent heat with "
               "respect to temperature "
            << approximated_dLw_dT
            << " and for computed derivative of water vapour latent heat "
               "with respect to temperature."
            << analytic_dLw_dT;
    }
}
