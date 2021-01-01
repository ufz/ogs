/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>
#include <sstream>

#include "TestMPL.h"
#include "Tests/TestTools.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/IdealGasLaw.h"
#include "MaterialLib/PhysicalConstant.h"

TEST(MaterialPropertyLib, IdealGasLawOfPurePhase)
{
    const double pressure_norm = 101325.0;    // Pa
    const double temperature_norm = 273.15;   // K
    const double molar_mass_air = 0.02905;  // kg/mol
    const auto R = MaterialLib::PhysicalConstant::IdealGasConstant;

    const double density_norm_air =
        pressure_norm * molar_mass_air / R / temperature_norm;

    const double d_rho_dT_air = -density_norm_air / temperature_norm;
    const double d_rho_dp_air = density_norm_air / pressure_norm;

    const double d_rho2_dp2_air = 0.;
    const double d_rho2_dT2_air =
        2. * density_norm_air / temperature_norm / temperature_norm;
    const double d_rho2_dTdp_air =
        -density_norm_air / temperature_norm / pressure_norm;

    std::stringstream m;
    m << "<medium>\n";
    m << "<phases><phase>\n";
    m << "  <type>Gas</type>\n";
    m << "  <properties>\n";
    m << "    <property>\n";
    m << "      <name>density</name>\n";
    m << "      <type>IdealGasLaw</type>\n";
    m << "    </property>\n";
    m << "    <property>\n";
    m << "      <name>molar_mass</name>\n";
    m << "      <type>Constant</type>\n";
    m << "      <value>" << molar_mass_air << "</value>\n";
    m << "    </property>\n";
    m << "  </properties>\n";
    m << "</phase></phases>\n";
    m << "<properties></properties>\n";
    m << "</medium>\n";

    auto const& medium = Tests::createTestMaterial(m.str());
    auto const& gas_phase = medium->phase("Gas");

    MaterialPropertyLib::VariableArray variable_array;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)] = pressure_norm;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = temperature_norm;

    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    auto const density =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(variable_array, pos, time, dt);

    auto const d_rho_dT =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(variable_array,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, time, dt);

    auto const d_rho_dp =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                variable_array, MaterialPropertyLib::Variable::phase_pressure,
                pos, time, dt);

    auto const d_rho2_dp2 =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template d2Value<double>(
                variable_array, MaterialPropertyLib::Variable::phase_pressure,
                MaterialPropertyLib::Variable::phase_pressure, pos, time, dt);

    auto const d_rho2_dT2 =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template d2Value<double>(
                variable_array, MaterialPropertyLib::Variable::temperature,
                MaterialPropertyLib::Variable::temperature, pos, time, dt);

    auto const d_rho2_dTdp =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template d2Value<double>(
                variable_array, MaterialPropertyLib::Variable::temperature,
                MaterialPropertyLib::Variable::phase_pressure, pos, time, dt);

    auto const d_rho2_dpdT =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template d2Value<double>(
                variable_array, MaterialPropertyLib::Variable::phase_pressure,
                MaterialPropertyLib::Variable::temperature, pos, time, dt);

    ASSERT_NEAR(density_norm_air, density, 1.e-10);
    ASSERT_NEAR(d_rho_dT_air, d_rho_dT, 1.e-10);
    ASSERT_NEAR(d_rho_dp_air, d_rho_dp, 1.e-10);
    ASSERT_NEAR(d_rho2_dp2_air, d_rho2_dp2, 1.e-10);
    ASSERT_NEAR(d_rho2_dT2_air, d_rho2_dT2, 1.e-10);
    ASSERT_NEAR(d_rho2_dTdp_air, d_rho2_dTdp, 1.e-10);
    ASSERT_EQ(d_rho2_dTdp, d_rho2_dpdT);
}
