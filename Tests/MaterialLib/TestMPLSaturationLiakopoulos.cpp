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
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationLiakopoulos.h"

TEST(MaterialPropertyLib, SaturationLiakopoulos)
{
    std::stringstream m;
    m << "<medium>\n";
    m << "<phases></phases>\n";
    m << "<properties>\n";
    m << "  <property>\n";
    m << "    <name>saturation</name>\n";
    m << "    <type>SaturationLiakopoulos</type>\n";
    m << "  </property> \n";
    m << "</properties>\n";
    m << "</medium>\n";

    auto const& medium = Tests::createTestMaterial(m.str());

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::array<double, 17> p_cap = {
        -1.000000000E+00, 1.000000000E+00, 2.000000000E+00, 4.000000000E+00,
        1.000000000E+01,  2.000000000E+01, 4.000000000E+01, 1.000000000E+02,
        2.000000000E+02,  4.000000000E+02, 1.000000000E+03, 2.000000000E+03,
        4.000000000E+03,  1.000000000E+04, 2.340000000E+04, 2.340265994E+04,
        2.340266000E+04};
    std::array<double, 17> s_L_ref = {
        1.0000000000E+00, 9.9999999998E-01, 9.9999999989E-01, 9.9999999943E-01,
        9.9999999472E-01, 9.9999997157E-01, 9.9999984703E-01, 9.9999858502E-01,
        9.9999238585E-01, 9.9995902751E-01, 9.9962098975E-01, 9.9796050953E-01,
        9.8902530632E-01, 8.9848015302E-01, 2.0022074567E-01, 2.0000000002E-01,
        2.0000000000E-01};

    std::array<double, 17> ds_L_dp_cap_ref = {
        0.0000000000E+00,  -4.7883043800E-11, -1.2883162371E-10,
        -3.4662765668E-10, -1.2825719853E-09, -3.4508213823E-09,
        -9.2846002786E-09, -3.4354351082E-08, -9.2432027708E-08,
        -2.4869279952E-07, -9.2019898433E-07, -2.4758394596E-06,
        -6.6613646984E-06, -2.4648003648E-05, -8.2982224427E-05,
        -8.2995693863E-05, -8.2995694166E-05};

    std::array<double, 17> d2s_L_dp_cap2_ref = {
        0.0000000000E+00,  -6.8372198242E-11, -9.1979337747E-11,
        -1.2373740774E-10, -1.8313845379E-10, -2.4637139259E-10,
        -3.3143701844E-10, -4.9054577909E-10, -6.5991846182E-10,
        -8.8777112110E-10, -1.3139521297E-09, -1.7676255822E-09,
        -2.3779406632E-09, -3.5194884408E-09, -5.0636888145E-09,
        -5.0639351070E-09, -5.0639351125E-09};

    for (std::size_t idx = 0; idx < p_cap.size(); idx++)
    {
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = p_cap[idx];

        auto s_L =
            medium->property(MaterialPropertyLib::PropertyType::saturation)
                .template value<double>(variable_array, pos, time, dt);
        auto ds_L_dp_cap =
            medium->property(MaterialPropertyLib::PropertyType::saturation)
                .template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::capillary_pressure, pos,
                    time, dt);
        auto d2s_L_dp_cap2 =
            medium->property(MaterialPropertyLib::PropertyType::saturation)
                .template d2Value<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::capillary_pressure,
                    MaterialPropertyLib::Variable::capillary_pressure, pos,
                    time, dt);

        ASSERT_NEAR(s_L, s_L_ref[idx], 1.e-10);
        ASSERT_NEAR(ds_L_dp_cap, ds_L_dp_cap_ref[idx], 1.e-10);
        ASSERT_NEAR(d2s_L_dp_cap2, d2s_L_dp_cap2_ref[idx], 1.e-10);
    }
}
