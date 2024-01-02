/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 20, 2022
 */

#include <gtest/gtest.h>

#include <sstream>
/* Keep for debugging
#include <fstream>
#include <iostream>
*/

#include "MaterialLib/MPL/Properties/TemperatureDependentFraction.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

struct ParameterSet
{
    double const k = 1;
    double const T_c = 273.15;  // K
};

std::unique_ptr<MaterialPropertyLib::Medium> createMyMedium()
{
    std::stringstream prj;
    ParameterSet water_ice;

    prj << "<medium>\n";
    prj << "    <properties>\n";
    prj << "         <property>\n";
    prj << "              <name>porosity</name>\n";
    prj << "              <type>Constant</type>\n";
    prj << "              <value>0.46</value>\n";
    prj << "         </property>\n";
    prj << "         <property>\n";
    prj << "              <name>volume_fraction</name>\n";
    prj << "              <type>TemperatureDependentFraction</type>\n";
    prj << "              <steepness>" << water_ice.k << "</steepness>\n";
    prj << "              <characteristic_temperature>" << water_ice.T_c
        << "</characteristic_temperature>\n";
    prj << "          </property>\n";
    prj << "    </properties>\n";
    prj << "</medium>\n";

    return Tests::createTestMaterial(prj.str());
}

/* Keep for debugging
TEST(MaterialPropertyLib, TemperatureDependentFraction_values_out)
{
    // arbitrary values for space and time
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    MPL::VariableArray vars;
    ParameterSet water_ice;

    double T_c = water_ice.T_c;
    double T_s = T_c-10;
    double T_e = T_c+10;
    size_t n_T = 400;
    double dT = (T_e - T_s) / n_T;

    auto const medium = createMyMedium();
    auto const& fraction =
        medium->property(MaterialPropertyLib::PropertyType::volume_fraction);

    std::vector<double> temperatures(n_T+1);
    std::generate(temperatures.begin(), temperatures.end(),
                  [T = T_s - dT, dT]() mutable { return T += dT; });

    std::ofstream fout("/home/cbs/data_frozen_fraction.txt");
    fout << "T" << '\t' << "f" << '\t' << "f,T" << '\t' << "f,TT" << '\n';
    for (double T : temperatures)
    {
        vars.temperature = T;

        auto const f = std::get<double>(fraction.value(vars, pos, time, dt));
        auto const fT = std::get<double>(
            fraction.dValue(vars, MPL::Variable::temperature, pos, time, dt));
        auto const fTT = std::get<double>(
            fraction.d2Value(vars, MPL::Variable::temperature,
                             MPL::Variable::temperature, pos, time, dt));
        if (!fout)
        {
            std::cout << T - T_c << '\t' << f << '\t' << fT << '\t' << fTT
                      << '\n';
        }
        else
        {
            fout << T-T_c << '\t' << f << '\t' << fT << '\t' << fTT << '\n';
        }
    }
    fout.close();
}
*/

TEST(MaterialPropertyLib, TemperatureDependentFraction_value_atTc)
{
    // arbitrary values for space and time
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    MPL::VariableArray vars;
    ParameterSet water_ice;
    // create medium with properties
    auto const medium = createMyMedium();

    vars.temperature = water_ice.T_c;

    auto const phi =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(vars, pos, time, dt);

    // calculate the temperature dependent (e.g. frozen) part of the pore space
    auto const pfr_expected = phi * 0.5;
    auto const pfr =
        medium->property(MaterialPropertyLib::PropertyType::volume_fraction)
            .template value<double>(vars, pos, time, dt);

    auto const relativeError = std::fabs((pfr_expected - pfr) / pfr_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected apparent heat capacity " << pfr_expected
        << " and for actual apparent heat capacity " << pfr;
}

TEST(MaterialPropertyLib, TemperatureDependentFraction_dValue_atTc)
{
    // arbitrary values for space and time
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    MPL::VariableArray vars;
    ParameterSet water_ice;
    // create medium with properties
    auto const medium = createMyMedium();

    vars.temperature = water_ice.T_c;

    auto const phi =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(vars, pos, time, dt);

    // calculate temperature-derivative of temperature dependent (e.g. frozen)
    // volume fraction
    auto const dpfr_dT_expected = -phi * water_ice.k * 0.25;
    auto const dpfr_dT =
        medium->property(MaterialPropertyLib::PropertyType::volume_fraction)
            .template dValue<double>(vars, MPL::Variable::temperature, pos,
                                     time, dt);

    auto const relativeError =
        std::fabs((dpfr_dT_expected - dpfr_dT) / dpfr_dT_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected apparent heat capacity " << dpfr_dT_expected
        << " and for actual apparent heat capacity " << dpfr_dT;
}
