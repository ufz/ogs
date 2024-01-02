/*!
   \file
   \brief Test class for water enthalpy in IAPWSIF97 region1 model.

   \author Chaofan Chen
   \date March 2023

   \copyright
    Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <array>
#include <memory>

#include "MaterialLib/MPL/Properties/Enthalpy/CreateWaterEnthalpyIAPWSIF97Region1.h"
#include "TestMPL.h"

TEST(Material, checkWaterEnthalpyIAPWSIF97Region1)
{
    const char xml[] =
        "<property>"
        "   <name>enthalpy</name>"
        "   <type>WaterEnthalpyIAPWSIF97Region1</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterEnthalpyIAPWSIF97Region1);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::array const pressures = {0., 3.e6, 80.e6, 3.e6};
    std::array const temperatures = {273.15, 300., 300., 500.};

    double const expected_h[] = {-0.042208551, 0.115325859e3, 0.184134183e3,
                                 0.975496445e3};

    for (std::size_t i = 0; i < pressures.size(); i++)
    {
        double const T_i = temperatures[i];
        double const p_i = pressures[i];

        variable_array.liquid_phase_pressure = p_i;
        variable_array.temperature = T_i;

        ASSERT_NEAR(
            expected_h[i],
            property.template value<double>(variable_array, pos, t, dt) / 1e3,
            1.e-6);

        double const dT = 1.0e-5;
        variable_array.temperature = T_i - dT;
        double const h0 =
            property.template value<double>(variable_array, pos, t, dt);

        variable_array.temperature = T_i + dT;
        double const h1 =
            property.template value<double>(variable_array, pos, t, dt);

        double const approximated_dhdT = 0.5 * (h1 - h0) / dT;
        double const analytic_dhdT = property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::temperature, pos, t,
            dt);

        ASSERT_LE(std::fabs(approximated_dhdT - analytic_dhdT), 2e-4)
            << "for expected derivative of water enthalpy in IAPWSIF97 region "
               "1 with "
               "respect to temperature "
            << approximated_dhdT;
    }
}
