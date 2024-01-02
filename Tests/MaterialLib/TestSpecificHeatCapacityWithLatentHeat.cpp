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

#include "TestMPL.h"
#include "Tests/TestTools.h"

struct IceWaterRockParameters
{
    double const k = 1;
    double const T_c = 273.15;   // K
    double const rho_I = 900;    // kg/m³
    double const rho_W = 1000;   // kg/m³
    double const rho_R = 3000;   // kg/m³
    double const kap_I = 2.37;   // W/m/K
    double const kap_W = 0.54;   // W/m/K
    double const kap_R = 3.00;   // W/m/K
    double const cp_I = 2052;    // J/kg/K
    double const cp_W = 4186;    // J/kg/K
    double const cp_R = 2500;    // J/kg/K
    double const L_IW = 334.e3;  // J/kg
};

std::unique_ptr<MaterialPropertyLib::Medium> createMyMedium(double L_IW,
                                                            double porosity)
{
    std::stringstream prj;
    IceWaterRockParameters water_ice_rock;

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
    prj << "           <name>specific_heat_capacity</name>\n";
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
    prj << "           <name>specific_heat_capacity</name>\n";
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
    prj << "           <name>specific_heat_capacity</name>\n";
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
    prj << "        <value>" << porosity << "</value>\n";
    prj << "      </property>\n";
    prj << "      <property>\n";
    prj << "        <name>volume_fraction</name>\n";
    prj << "        <type>TemperatureDependentFraction</type>\n";
    prj << "        <steepness>" << water_ice_rock.k << "</steepness>\n";
    prj << "        <characteristic_temperature>" << water_ice_rock.T_c
        << "</characteristic_temperature>\n";
    prj << "      </property>\n";
    prj << "       <property>\n";
    prj << "         <name>density</name>\n";
    prj << "         <type>VolumeFractionAverage</type>\n";
    prj << "       </property>\n";
    prj << "      <property>\n";
    prj << "        <name>thermal_conductivity</name>\n";
    prj << "        <type>VolumeFractionAverage</type>\n";
    prj << "      </property>\n";
    prj << "      <property>\n";
    prj << "        <name>specific_heat_capacity</name>\n";
    prj << "        <type>SpecificHeatCapacityWithLatentHeat</type>\n";
    prj << "        <specific_latent_heat>" << L_IW
        << "</specific_latent_heat>\n";
    prj << "      </property>\n";
    prj << "    </properties>\n";
    prj << "</medium>\n";

    return Tests::createTestMaterial(prj.str());
}

// Test for trivial case of zero latent heat at the critical temperature
TEST(MaterialPropertyLib, SpecificHeatCapacityWithLatentHeat_trivial)
{
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    IceWaterRockParameters water_ice_rock;

    auto const zero_latent_heat_IW = 0.0;
    auto const full_porosity = 1.0;
    auto const& medium = createMyMedium(zero_latent_heat_IW, full_porosity);

    // at the critical temperature the frozen and the liquid phase have
    // the same volume fractions, there is no solid involved
    auto const rho_mix = 0.5 * (water_ice_rock.rho_I + water_ice_rock.rho_W);
    auto const Cvol_mix = 0.5 * (water_ice_rock.rho_I * water_ice_rock.cp_I +
                                 water_ice_rock.rho_W * water_ice_rock.cp_W);

    vars.temperature = water_ice_rock.T_c;
    vars.density = rho_mix;

    auto const Capp =
        medium
            ->property(
                MaterialPropertyLib::PropertyType::specific_heat_capacity)
            .template value<double>(vars, pos, t, dt);
    auto const Capp_expected = Cvol_mix / rho_mix;
    auto const relativeError =
        std::fabs((Capp_expected - Capp) / Capp_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected apparent heat capacity " << Capp_expected
        << " and for actual apparent heat capacity " << Capp;
}

// Test for an arbitrary set of values at the critical temperature
TEST(MaterialPropertyLib, SpecificHeatCapacityWithLatentHeat_atTc)
{
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    IceWaterRockParameters water_ice_rock;

    auto const phi = 0.46;
    auto const& medium = createMyMedium(water_ice_rock.L_IW, phi);

    // at the critical temperature the frozen and the liquid phase have same
    // volume fractions of 0.5
    auto const rho_mix =
        phi * 0.5 * (water_ice_rock.rho_I + water_ice_rock.rho_W) +
        (1 - phi) * water_ice_rock.rho_R;

    auto const Cvol_mix =
        phi * 0.5 *
            (water_ice_rock.rho_I * water_ice_rock.cp_I +
             water_ice_rock.rho_W * water_ice_rock.cp_W) +
        (1 - phi) * water_ice_rock.rho_R * water_ice_rock.cp_R;

    vars.temperature = water_ice_rock.T_c;
    vars.density = rho_mix;

    auto const Capp =
        medium->property(MPL::PropertyType::specific_heat_capacity)
            .template value<double>(vars, pos, t, dt);
    auto const Capp_expected =
        (Cvol_mix + water_ice_rock.rho_I * water_ice_rock.L_IW * phi *
                        water_ice_rock.k / 4) /
        rho_mix;
    auto const relativeError =
        std::fabs((Capp_expected - Capp) / Capp_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected apparent heat capacity " << Capp_expected
        << " and for actual apparent heat capacity " << Capp;
}

// Test for an arbitrary set of values far below the critical temperature
TEST(MaterialPropertyLib, SpecificHeatCapacityWithLatentHeat_belowTc)
{
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    IceWaterRockParameters water_ice_rock;

    auto const full_porosity = 1.0;
    auto const& medium = createMyMedium(water_ice_rock.L_IW, full_porosity);

    vars.temperature = 0.5 * water_ice_rock.T_c;
    vars.density = water_ice_rock.rho_I;

    auto const Capp =
        medium->property(MPL::PropertyType::specific_heat_capacity)
            .template value<double>(vars, pos, t, dt);
    auto const Capp_expected = water_ice_rock.cp_I;
    auto const relativeError =
        std::fabs((Capp_expected - Capp) / Capp_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected apparent heat capacity " << Capp_expected
        << " and for actual apparent heat capacity " << Capp;
}

// Test for an arbitrary set of values far above the critical temperature
TEST(MaterialPropertyLib, SpecificHeatCapacityWithLatentHeat_aboveTc)
{
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    IceWaterRockParameters water_ice_rock;

    auto const full_porosity = 1.0;
    auto const& medium = createMyMedium(water_ice_rock.L_IW, full_porosity);

    vars.temperature = 1.25 * water_ice_rock.T_c;
    vars.density = water_ice_rock.rho_W;

    auto const Capp =
        medium->property(MPL::PropertyType::specific_heat_capacity)
            .template value<double>(vars, pos, t, dt);
    auto const Capp_expected = water_ice_rock.cp_W;
    auto const relativeError =
        std::fabs((Capp_expected - Capp) / Capp_expected);

    ASSERT_LE(relativeError, 1e-10)
        << "for expected apparent heat capacity " << Capp_expected
        << " and for actual apparent heat capacity " << Capp;
}
