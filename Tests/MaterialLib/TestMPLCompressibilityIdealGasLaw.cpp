/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MaterialLib/MPL/Properties/CompressibilityIdealGasLaw.h"

TEST(MaterialPropertyLib, CompressibilityIdealGasLaw)
{
    const double pressure_norm = 101325.0;    // Pa

    const double compressibility_norm_air = 1.0 / pressure_norm;

    const double d_beta_dp_air = compressibility_norm_air / pressure_norm;

    const double d_beta_dp2_air = compressibility_norm_air /
        (pressure_norm * pressure_norm);

    std::stringstream m;
    m << "<medium>\n";
    m << "<phases><phase>\n";
    m << "  <type>Gas</type>\n";
    m << "  <properties>\n";
    m << "    <property>\n";
    m << "      <name>compressibility</name>\n";
    m << "      <type>CompressibilityIdealGasLaw</type>\n";
    m << "    </property>\n";
    m << "  </properties>\n";
    m << "</phase></phases>\n";
    m << "<properties></properties>\n";
    m << "</medium>\n";

    auto const& medium = createTestMaterial(m.str());
    auto const& gas_phase = medium->phase("Gas");

    MaterialPropertyLib::VariableArray variable_array;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)] = pressure_norm;

    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();

    auto const beta =
        gas_phase.property(MaterialPropertyLib::PropertyType::compressibility)
            .template value<double>(variable_array, pos, time);

    auto const d_beta_dp =
        gas_phase.property(MaterialPropertyLib::PropertyType::compressibility)
            .template dValue<double>(
                variable_array, MaterialPropertyLib::Variable::phase_pressure,
                pos, time);

    auto const d_beta_dp2 =
        gas_phase.property(MaterialPropertyLib::PropertyType::density)
            .template d2Value<double>(
                variable_array, MaterialPropertyLib::Variable::phase_pressure,
                MaterialPropertyLib::Variable::phase_pressure, pos, time);

    ASSERT_NEAR(compressibility_norm_air, beta, 1.e-10);
    ASSERT_NEAR(d_beta_dp_air, d_beta_dp, 1.e-10);
    ASSERT_NEAR(d_beta_dp2_air, d_beta_dp2, 1.e-10);
}
