/**
 *  \file
 *  \brief Test thermal conductivity of water IAPWS97 model
 *
 *  \copyright
 *   Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/ThermalConductivity/CreateWaterThermalConductivityIAPWS.h"
#include "TestMPL.h"

TEST(Material, checkWaterThermalConductivityIAPWS_)
{
    const char xml[] =
        "<property>"
        "  <name>thermal_conductivity</name>"
        "  <type>WaterThermalConductivityIAPWS</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterThermalConductivityIAPWS);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // Test data provided on http://www.iapws.org/relguide/ThCond.pdf Table 4
    const double T[] = {298.15, 298.15, 298.15, 873.15};
    const double rho[] = {0, 998, 1200, 0};

    const double expected_lambda[] = {18.4341883e-3, 607.712868e-3,
                                      799.038144e-3, 79.1034659e-3};

    const double perturbation = 1.e-8;
    for (int i = 0; i < 4; i++)
    {
        // Test lambda
        variable_array.temperature = T[i];
        variable_array.density = rho[i];
        const double lambda =
            property.template value<double>(variable_array, pos, t, dt);
        ASSERT_NEAR(expected_lambda[i], lambda, 1.e-9);

        const double dlambda_dT = property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::temperature, pos, t,
            dt);
        const double dlambda_drho = property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::density, pos, t, dt);

        // Test dlambda/dT
        variable_array.temperature = T[i] + perturbation;
        double lambda1 =
            property.template value<double>(variable_array, pos, t, dt);
        ASSERT_NEAR((lambda1 - lambda) / perturbation, dlambda_dT, 8.e-6);

        // Test dlambda/drho
        variable_array.temperature = T[i];
        variable_array.density = rho[i] + perturbation;
        lambda1 = property.template value<double>(variable_array, pos, t, dt);

        ASSERT_NEAR((lambda1 - lambda) / perturbation, dlambda_drho, 8.e-6);
    }
}
