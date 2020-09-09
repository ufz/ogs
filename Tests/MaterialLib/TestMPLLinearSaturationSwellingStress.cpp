/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 14, 2020, 9:06 AM
 */

#include <gtest/gtest.h>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Properties/SwellingStress/CreateLinearSaturationSwellingStress.h"
#include "MaterialLib/MPL/Properties/SwellingStress/LinearSaturationSwellingStress.h"
#include "Tests/TestTools.h"

std::unique_ptr<MaterialPropertyLib::Property>
createLinearSaturationSwellingStressModel(const char xml[])
{
    auto const ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("property");
    return MaterialPropertyLib::createLinearSaturationSwellingStress(
        sub_config);
}

TEST(MaterialPropertyLib, LinearSaturationSwellingStress)
{
    const char xml[] =
        "<property>"
        "    <name>swelling_stress_rate</name>"
        "    <type>LinearSaturationSwellingStress</type>"
        "    <coefficient>1.e+6 </coefficient>"
        "    <reference_saturation> 0.65 </reference_saturation>"
        "</property>";

    auto const swelling_stress_rate =
        createLinearSaturationSwellingStressModel(xml);

    double const coefficient = 1.0e+6;
    double const dt = 0.1;

    // Saturation is larger than the reference saturation saturation.
    double const S0 = 0.65;
    double const S1 = 0.7;
    const double dS = S1 - S0;

    MaterialPropertyLib::VariableArray variable_array;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation_rate)] = dS / dt;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = S1;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double d_sw_expected = coefficient * dS;
    double d_sw = std::get<double>(
        swelling_stress_rate->value(variable_array, pos, time, dt));
    ASSERT_LE(std::fabs(d_sw_expected - d_sw), 1e-19)
        << "for expected swelling stress rate" << d_sw_expected
        << " for computed swelling stress rate." << d_sw_expected;

    double dsw_dS = std::get<double>(swelling_stress_rate->dValue(
        variable_array, MaterialPropertyLib::Variable::liquid_saturation, pos,
        time, dt));
    double dsw_dS_expected = coefficient;
    ASSERT_LE(std::fabs(dsw_dS_expected - dsw_dS), 1e-19)
        << "for expected dsw/dS" << d_sw_expected << " for computed dsw/dS."
        << d_sw_expected;

    // Saturation is smaller than the reference saturation saturation.
    double const S2 = 0.3;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = S2;

    d_sw = std::get<double>(
        swelling_stress_rate->value(variable_array, pos, time, dt));
    d_sw_expected = 0.0;
    ASSERT_LE(std::fabs(d_sw_expected - d_sw), 1e-19)
        << "for expected swelling stress rate" << d_sw_expected
        << " for computed swelling stress rate." << d_sw_expected;

    dsw_dS = std::get<double>(swelling_stress_rate->dValue(
        variable_array, MaterialPropertyLib::Variable::liquid_saturation, pos,
        time, dt));
    dsw_dS_expected = 0.0;
    ASSERT_LE(std::fabs(dsw_dS_expected - dsw_dS), 1e-19)
        << "for expected dsw/dS" << d_sw_expected << " for computed dsw/dS."
        << d_sw_expected;
}
