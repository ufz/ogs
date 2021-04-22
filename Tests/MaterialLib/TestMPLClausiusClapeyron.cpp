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
#include <string>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/ClausiusClapeyron.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

std::stringstream config(const double molar_mass,
                         const double T_triple,
                         const double T_critical,
                         const double T_ref,
                         const double p_triple,
                         const double p_critical,
                         const double p_ref,
                         const double dMdT,
                         const double dMdp)
{
    std::stringstream m;
    m << std::setprecision(16);
    m << "<medium>\n";
    m << "    <phases>\n";
    m << "        <phase>\n";
    m << "            <type>Gas</type>\n";
    m << "            <components>\n";
    m << "                <component>\n";
    m << "                    <name>W</name>\n";
    m << "                    <properties>\n";
    m << " <property>\n";
    m << "   <name>molar_mass</name>\n";
    m << "     <type>Linear</type>\n";
    m << "       <reference_value>" << molar_mass << " </reference_value>\n";
    m << "         <independent_variable>\n";
    m << "           <variable_name>temperature</variable_name>\n";
    m << "           <reference_condition>0</reference_condition>\n";
    m << "           <slope>" << dMdT << "</slope>\n";
    m << "         </independent_variable>\n";
    m << "         <independent_variable>\n";
    m << "           <variable_name>phase_pressure</variable_name>\n";
    m << "           <reference_condition>0.0</reference_condition>\n";
    m << "           <slope>" << dMdp << "</slope>\n";
    m << "         </independent_variable>\n";
    m << " </property>\n";
    m << "                        <property>\n";
    m << "                            <name>vapour_pressure</name>\n";
    m << "                            <type>ClausiusClapeyron</type>\n";
    m << "<triple_temperature>" << T_triple << "</triple_temperature>\n";
    m << "<critical_temperature>" << T_critical << "</critical_temperature>\n";
    m << "<reference_temperature>" << T_ref << "</reference_temperature>\n";
    m << "<triple_pressure>" << p_triple << "</triple_pressure>\n";
    m << "<critical_pressure>" << p_critical << "</critical_pressure>\n";
    m << "<reference_pressure>" << p_ref << "</reference_pressure>\n";
    m << "                        </property>\n";
    m << "                    </properties>\n";
    m << "                </component>\n";
    m << "            </components>\n";
    m << "            <properties>\n";
    m << "            </properties>\n";
    m << "        </phase>\n";
    m << "    </phases>\n";
    m << "    <properties>\n";
    m << "    </properties>\n";
    m << "</medium>\n";
    return m;
}

const double molar_mass = 0.018015;
const double evaporation_enthalpy = 2.257e6;
const double T_triple = 273.16;
const double T_critical = 647.1;
const double T_ref = 373.15;
const double p_triple = 836.21849;
const double p_critical = 26016399.68391;
const double p_ref = 101325;
const double dMdT = 1.e-5;
const double dMdp = -1.e-8;

std::array<double, 29> const ref_temperature = {
    2.70000E+02, 2.72000E+02, 2.73000E+02, 2.73150E+02, 2.73160E+02,
    2.75000E+02, 2.80000E+02, 3.00000E+02, 3.20000E+02, 3.40000E+02,
    3.60000E+02, 3.73150E+02, 4.00000E+02, 4.50000E+02, 5.00000E+02,
    5.50000E+02, 6.00000E+02, 6.47100E+02, 6.50000E+02, 2.70000E+02,
    3.00000E+02, 3.30000E+02, 3.60000E+02, 3.90000E+02, 4.20000E+02,
    2.70000E+02, 3.00000E+02, 3.30000E+02, 3.60000E+02};

std::array<double, 29> const ref_liquid_pressure = {
    1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05,
    1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05,
    1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05,
    1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05, 1.01325E+05,
    1.50000E+05, 2.00000E+05, 5.00000E+05, 1.00000E+06, 5.00000E+06,
    1.01325E+05, 1.50000E+05, 2.00000E+05, 5.00000E+05};

std::array<double, 29> const ref_vapour_pressure = {
    8.36218490000000E+02, 8.36218490000000E+02, 8.36218490000000E+02,
    8.36218490000000E+02, 8.36218490000000E+02, 9.35003929341767E+02,
    1.28489641893467E+03, 4.12252674333030E+03, 1.14370551053552E+04,
    2.81490332971881E+04, 6.27015396656223E+04, 1.01325000000000E+05,
    2.44852602635238E+05, 9.57415384221067E+05, 2.85381686957823E+06,
    6.98277711336414E+06, 1.47345726218355E+07, 2.60163996839100E+07,
    2.60163996839100E+07, 8.36218490000000E+02, 4.12894398711138E+03,
    1.82193325072170E+04, 6.28213193519972E+04, 1.77877844206557E+05,
    4.08809303341452E+05, 8.36218490000000E+02, 4.12894398711138E+03,
    1.82193325072170E+04, 6.28213193519972E+04};

std::array<double, 29> const ref_d_vapour_pressure_dT = {
    0.00000000000000E+00, 0.00000000000000E+00, 0.00000000000000E+00,
    0.00000000000000E+00, 5.48586853887215E+01, 6.05228838489413E+01,
    8.02335271069103E+01, 2.24315684282508E+02, 5.47138176563398E+02,
    1.19327831787700E+03, 2.37176618859269E+03, 3.56829135217340E+03,
    7.50820659017266E+03, 2.32230679272202E+04, 5.61410287431912E+04,
    1.13685276471403E+05, 2.01883331325872E+05, 3.06935462457948E+05,
    0.00000000000000E+00, 0.00000000000000E+00, 2.24555657648407E+02,
    8.18907637521856E+02, 2.36684655384022E+03, 5.68517643381087E+03,
    1.08201481763584E+04, 0.00000000000000E+00, 2.24555657648407E+02,
    8.18907637521856E+02, 2.36684655384022E+03};

std::array<double, 29> const ref_d_vapour_pressure_dp = {
    0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
    0.00000000000000E+00,  4.01150605544652E-05,  4.37340067538689E-05,
    5.60197502360705E-05,  1.31736092964623E-04,  2.48952060831370E-04,
    3.59680026912301E-04,  3.00157831341098E-04,  -0.00000000000000E+00,
    -2.15395807049717E-03, -2.14279182737438E-02, -9.48843943114251E-02,
    -2.94251561405101E-01, -7.30084785980384E-01, -1.44342865427858E+00,
    0.00000000000000E+00,  0.00000000000000E+00,  1.31941157158491E-04,
    3.12210527929795E-04,  3.00731227322962E-04,  -1.00717602020783E-03,
    -5.97625743435430E-03, 0.00000000000000E+00,  1.31941157158491E-04,
    3.12210527929795E-04,  3.00731227322962E-04};

TEST(MaterialPropertyLib, ClausiusClapeyron)
{
    std::stringstream m;

    for (std::size_t i = 0; i < ref_temperature.size(); i++)
    {
        MaterialPropertyLib::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const t = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();

        std::stringstream m = config(molar_mass, T_triple, T_critical, T_ref,
                                     p_triple, p_critical, p_ref, dMdT, dMdp);

        auto const& medium = Tests::createTestMaterial(m.str());
        auto const& gas_phase = medium->phase("Gas");

        const double p_LR = ref_liquid_pressure[i];
        const double T = ref_temperature[i];

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::phase_pressure)] = p_LR;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::enthalpy_of_evaporation)] =
            evaporation_enthalpy;

        auto const p_vap =
            gas_phase.component(0)
                .property(MaterialPropertyLib::PropertyType::vapour_pressure)
                .template value<double>(variable_array, pos, t, dt);

        auto const d_p_vap_dT =
            gas_phase.component(0)
                .property(MaterialPropertyLib::PropertyType::vapour_pressure)
                .template dValue<double>(
                    variable_array, MaterialPropertyLib::Variable::temperature,
                    pos, t, dt);

        auto const d_p_vap_dp =
            gas_phase.component(0)
                .property(MaterialPropertyLib::PropertyType::vapour_pressure)
                .template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::phase_pressure, pos, t, dt);

        ASSERT_NEAR(p_vap, ref_vapour_pressure[i], 1.e-7);
        ASSERT_NEAR(d_p_vap_dT, ref_d_vapour_pressure_dT[i], 1.e-7);
        ASSERT_NEAR(d_p_vap_dp, ref_d_vapour_pressure_dp[i], 1.e-7);
    }
}
