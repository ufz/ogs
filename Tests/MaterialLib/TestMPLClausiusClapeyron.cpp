/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/ClausiusClapeyron.h"
#include "MaterialLib/PhysicalConstant.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

static const double triple_temperature = 273.16;
static const double triple_pressure = 611.66;
static const double critical_temperature = 6.47114E+02;
static const double critical_pressure = 2.20640E+07;
static const double ref_temperature = 373.15;
static const double ref_pressure = 101325;

static const double molar_mass = 0.0180156;
static const double vapourisation_enthalpy = 2258000.;

std::string mediumDefinition()
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
    m << Tests::makeConstantPropertyElement("molar_mass", molar_mass);
    m << "                        <property>\n";
    m << "                            <name>vapour_pressure</name>\n";
    m << "                            <type>ClausiusClapeyron</type>\n";
    m << "<triple_temperature>" << triple_temperature
      << "</triple_temperature>\n";
    m << "<critical_temperature>" << critical_temperature
      << "</critical_temperature>\n";
    m << "<reference_temperature>" << ref_temperature
      << "</reference_temperature>\n";
    m << "<triple_pressure>" << triple_pressure << "</triple_pressure>\n";
    m << "<critical_pressure>" << critical_pressure << "</critical_pressure>\n";
    m << "<reference_pressure>" << ref_pressure << "</reference_pressure>\n";
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
    return m.str();
}

TEST(MaterialPropertyLib, ClausiusClapeyron)
{
    auto const& medium = Tests::createTestMaterial(mediumDefinition());
    auto const& gas_phase = medium->phase("Gas");
    auto const& vapour_component = gas_phase.component(0);

    MPL::VariableArray vars;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    vars[static_cast<int>(MPL::Variable::enthalpy_of_evaporation)] =
        vapourisation_enthalpy;

    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double T_min = triple_temperature - 20;
    const double T_max = critical_temperature + 20;

    for (double T = T_min; T <= T_max; T += 0.12345)
    {
        vars[static_cast<int>(MPL::Variable::temperature)] = T;

        auto const pVap =
            vapour_component.property(MPL::PropertyType::vapour_pressure)
                .template value<double>(vars, pos, time, dt);

        auto const dpVap_dT =
            vapour_component.property(MPL::PropertyType::vapour_pressure)
                .template dValue<double>(
                    vars, MPL::Variable::temperature, pos, time, dt);

        auto const value =
            ref_pressure * std::exp((1. / ref_temperature - 1. / T) *
                                    molar_mass * vapourisation_enthalpy / R);

        auto const pVap_ = T <= triple_temperature     ? triple_pressure
                           : T >= critical_temperature ? critical_pressure
                                                       : value;
        ASSERT_NEAR(pVap, pVap_, 1.e-10);
        // perturb temperature to check derivatives
        auto const eps = 5.e-4;
        vars[static_cast<int>(MPL::Variable::temperature)] = T + eps;

        auto const pVap_plus =
            vapour_component.property(MPL::PropertyType::vapour_pressure)
                .template value<double>(vars, pos, time, dt);

        vars[static_cast<int>(MPL::Variable::temperature)] = T - eps;

        auto const pVap_minus =
            vapour_component.property(MPL::PropertyType::vapour_pressure)
                .template value<double>(vars, pos, time, dt);

        auto const dpVap_dT_ = (pVap_plus - pVap_minus) / (2 * eps);

        ASSERT_NEAR(dpVap_dT, dpVap_dT_, 1.e-04);
    }
}
