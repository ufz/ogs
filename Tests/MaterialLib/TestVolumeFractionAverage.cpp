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

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "TestMPL.h"

struct ParameterSetIceWaterRock
{
    double const k = 1;
    double const T_c = 273.15;  // K
    double const rho_I = 900;   // kg/m³
    double const rho_W = 1000;  // kg/m³
    double const rho_R = 3000;  // kg/m³
    double const kap_I = 2.37;  // W/m/K
    double const kap_W = 0.54;  // W/m/K
    double const kap_R = 3.00;  // W/m/K
    double const cp_I = 2052;   // J/kg/K
    double const cp_W = 4186;   // J/kg/K
    double const cp_R = 2500;   // J/kg/K
};

TEST(MaterialPropertyLib, VolumeFractionAverage_Density)
{
    std::stringstream prj;
    ParameterSetIceWaterRock water_ice_rock;

    prj << "<medium>\n";
    prj << "  <phases>\n";
    prj << "    <phase>\n";
    prj << "      <type>AqueousLiquid</type>\n";
    prj << "      <properties>\n";
    prj << "         <property>\n";
    prj << "           <name>density</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.rho_W << "</value>\n";
    prj << "         </property>\n";
    prj << "      </properties>\n";
    prj << "    </phase>\n";
    prj << "    <phase>\n";
    prj << "      <type>FrozenLiquid</type>\n";
    prj << "      <properties>\n";
    prj << "         <property>\n";
    prj << "           <name>density</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.rho_I << "</value>\n";
    prj << "         </property>\n";
    prj << "      </properties>\n";
    prj << "    </phase>\n";
    prj << "    <phase>\n";
    prj << "      <type>Solid</type>\n";
    prj << "      <properties>\n";
    prj << "         <property>\n";
    prj << "           <name>density</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.rho_R << "</value>\n";
    prj << "         </property>\n";
    prj << "      </properties>\n";
    prj << "    </phase>\n";
    prj << "  </phases>\n";
    prj << "    <properties>\n";
    prj << "       <property>\n";
    prj << "         <name>porosity</name>\n";
    prj << "         <type>Constant</type>\n";
    prj << "         <value>0.33</value>\n";
    prj << "       </property>\n";
    prj << "       <property>\n";
    prj << "         <name>volume_fraction</name>\n";
    prj << "         <type>TemperatureDependentFraction</type>\n";
    prj << "         <steepness>" << water_ice_rock.k << "</steepness>\n";
    prj << "         <characteristic_temperature>" << water_ice_rock.T_c
        << "</characteristic_temperature>\n";
    prj << "       </property>\n";
    prj << "       <property>\n";
    prj << "         <name>density</name>\n";
    prj << "         <type>VolumeFractionAverage</type>\n";
    prj << "       </property>\n";
    prj << "    </properties>\n";
    prj << "</medium>\n";

    auto const medium = Tests::createTestMaterial(prj.str());

    // ignored arguments: value must not have any effect -> NaN
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    // relevant argument:
    MaterialPropertyLib::VariableArray vars;
    vars.temperature = water_ice_rock.T_c;

    auto const& phi =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(vars, pos, time, dt);

    auto rho_m = medium->property(MaterialPropertyLib::PropertyType::density)
                     .template value<double>(vars, pos, time, dt);

    auto const rho_expected =
        (1 - phi) * water_ice_rock.rho_R +
        phi * 0.5 * (water_ice_rock.rho_W + water_ice_rock.rho_I);
    auto const relativeError = std::fabs((rho_expected - rho_m) / rho_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected mixture density " << rho_expected
        << " and for actual mixture density " << rho_m;
}

TEST(MaterialPropertyLib, VolumeFractionAverage_ThermalConductivity)
{
    std::stringstream prj;
    ParameterSetIceWaterRock water_ice_rock;

    prj << "<medium>\n";
    prj << "  <phases>\n";
    prj << "    <phase>\n";
    prj << "      <type>AqueousLiquid</type>\n";
    prj << "      <properties>\n";
    prj << "         <property>\n";
    prj << "           <name>density</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.rho_W << "</value>\n";
    prj << "         </property>\n";
    prj << "         <property>\n";
    prj << "           <name>heat_capacity</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.cp_W << "</value>\n";
    prj << "         </property>\n";
    prj << "         <property>\n";
    prj << "           <name>thermal_conductivity</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.kap_W << "</value>\n";
    prj << "         </property>\n";
    prj << "      </properties>\n";
    prj << "    </phase>\n";
    prj << "    <phase>\n";
    prj << "      <type>FrozenLiquid</type>\n";
    prj << "      <properties>\n";
    prj << "         <property>\n";
    prj << "           <name>density</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.rho_I << "</value>\n";
    prj << "         </property>\n";
    prj << "         <property>\n";
    prj << "           <name>heat_capacity</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.cp_I << "</value>\n";
    prj << "         </property>\n";
    prj << "         <property>\n";
    prj << "           <name>thermal_conductivity</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.kap_I << "</value>\n";
    prj << "         </property>\n";
    prj << "      </properties>\n";
    prj << "    </phase>\n";
    prj << "    <phase>\n";
    prj << "      <type>Solid</type>\n";
    prj << "      <properties>\n";
    prj << "         <property>\n";
    prj << "           <name>density</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.rho_R << "</value>\n";
    prj << "         </property>\n";
    prj << "         <property>\n";
    prj << "           <name>heat_capacity</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.cp_R << "</value>\n";
    prj << "         </property>\n";
    prj << "         <property>\n";
    prj << "           <name>thermal_conductivity</name>\n";
    prj << "           <type>Constant</type>\n";
    prj << "           <value>" << water_ice_rock.kap_R << "</value>\n";
    prj << "         </property>\n";
    prj << "      </properties>\n";
    prj << "    </phase>\n";
    prj << "  </phases>\n";
    prj << "    <properties>\n";
    prj << "      <property>\n";
    prj << "        <name>porosity</name>\n";
    prj << "        <type>Constant</type>\n";
    prj << "        <value>1.0</value>\n";
    prj << "      </property>\n";
    prj << "      <property>\n";
    prj << "        <name>volume_fraction</name>\n";
    prj << "        <type>TemperatureDependentFraction</type>\n";
    prj << "        <steepness>" << water_ice_rock.k << "</steepness>\n";
    prj << "        <characteristic_temperature>" << water_ice_rock.T_c
        << "</characteristic_temperature>\n";
    prj << "      </property>\n";
    prj << "      <property>\n";
    prj << "        <name>thermal_conductivity</name>\n";
    prj << "        <type>VolumeFractionAverage</type>\n";
    prj << "      </property>\n";
    prj << "    </properties>\n";
    prj << "</medium>\n";

    auto const medium = Tests::createTestMaterial(prj.str());

    // ignored arguments: value must not have any effect -> NaN
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    // relevant argument:
    MaterialPropertyLib::VariableArray vars;
    vars.temperature = water_ice_rock.T_c;

    auto const& phi =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(vars, pos, time, dt);

    auto kap_m =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(vars, pos, time, dt);

    auto const kap_expected =
        (1 - phi) * water_ice_rock.kap_R +
        phi * 0.5 * (water_ice_rock.kap_W + water_ice_rock.kap_I);

    auto const relativeError = std::fabs((kap_expected - kap_m) / kap_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected mixture thermal conductivity " << kap_expected
        << " and for actual mixture thermal conductivity " << kap_m;
}