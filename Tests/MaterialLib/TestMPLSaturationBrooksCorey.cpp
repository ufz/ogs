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
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationBrooksCorey.h"

TEST(MaterialPropertyLib, SaturationBrooksCorey)
{
    const double ref_lambda = 2.5;
    const double ref_residual_liquid_saturation = 0.12;
    const double ref_residual_gas_saturation = 0.06;
    const double ref_entry_pressure = 5678.54;

    const double max_saturation = 1. - ref_residual_gas_saturation;

    std::stringstream m;
    m << "<medium>\n";
    m << "<phases></phases>\n";
    m << "<properties>\n";
    m << "  <property>\n";
    m << "    <name>saturation</name>\n";
    m << "    <type>SaturationBrooksCorey</type>\n";
    m << "    <residual_liquid_saturation>" << ref_residual_liquid_saturation
      << "</residual_liquid_saturation>\n";
    m << "    <residual_gas_saturation>" << ref_residual_gas_saturation
      << "</residual_gas_saturation>\n";
    m << "    <lambda>" << ref_lambda << "</lambda>\n";
    m << "    <entry_pressure>" << ref_entry_pressure << "</entry_pressure>\n";
    m << "  </property> \n";
    m << "</properties>\n";
    m << "</medium>\n";

    auto const& medium = Tests::createTestMaterial(m.str());

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    for (double p_cap = 1.0; p_cap < 1.0e10; p_cap *= 1.5)
    {
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = p_cap;

        auto s_L =
            medium->property(MaterialPropertyLib::PropertyType::saturation)
                .template value<double>(variable_array, pos, t, dt);
        auto ds_L_dp_cap =
            medium->property(MaterialPropertyLib::PropertyType::saturation)
                .template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::capillary_pressure, pos, t,
                    dt);

        const double s_eff =
            std::pow(ref_entry_pressure / std::max(p_cap, ref_entry_pressure),
                     ref_lambda);
        const double s_ref =
            s_eff * (max_saturation - ref_residual_liquid_saturation) +
            ref_residual_liquid_saturation;
        const double ds_eff_dpc =
            -ref_lambda / std::max(p_cap, ref_entry_pressure) * s_ref;
        const double ds_L_ds_eff =
            1. / (max_saturation - ref_residual_liquid_saturation);
        const double ds_L_dpc = ds_L_ds_eff * ds_eff_dpc;

        ASSERT_NEAR(s_L, s_ref, 1.e-10);
        ASSERT_NEAR(ds_L_dp_cap, ds_L_dpc, 1.e-10);
    }
}
