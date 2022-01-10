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
#include "MaterialLib/MPL/Properties/EffectiveThermalConductivityPorosityMixing.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "TestMPL.h"

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingLiquidSolid1D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>AqueousLiquid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.923</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 1);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 1.0;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<1>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 1, 1.e-10);
}
TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingLiquidSolid2D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>AqueousLiquid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.923</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 2);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 1.0;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<2>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 2, 1.e-10);
}
TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingLiquidSolid3D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>AqueousLiquid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.923</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 3);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 1.0;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 3, 1.e-10);
}

TEST(MaterialPropertyLib, EffectiveThermalConductivityPorosityMixingGasSolid1D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>Gas</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.923</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 1);
    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 0.0;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<1>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 1, 1.e-10);
}

TEST(MaterialPropertyLib, EffectiveThermalConductivityPorosityMixingGasSolid2D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>Gas</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.923</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 2);
    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 0.0;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<2>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 2, 1.e-10);
}

TEST(MaterialPropertyLib, EffectiveThermalConductivityPorosityMixingGasSolid3D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>Gas</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.923</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 3);
    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 0.0;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 3, 1.e-10);
}

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingGasLiquidSolid1D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>Gas</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>AqueousLiquid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.456</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.789</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 1);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 0.3;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<1>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.721068 * 1, 1.e-10);
}

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingGasLiquidSolid2D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>Gas</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>AqueousLiquid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.456</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.789</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 2);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 0.3;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<2>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.721068 * 2, 1.e-10);
}

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingGasLiquidSolid3D)
{
    std::string m =
        "<medium>"
        "<phases><phase><type>Gas</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.123</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>AqueousLiquid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.456</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.789</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium = Tests::createTestMaterial(std::move(m), 3);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 0.3;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.721068 * 3, 1.e-10);
}

TEST(MaterialPropertyLib, EffectiveThermalConductivityPorosityMixingRot90deg)
{
    ParameterLib::ConstantParameter<double> const a{"a", {1., 0., 0.}};
    ParameterLib::ConstantParameter<double> const b{"a", {0., 1., 0.}};
    ParameterLib::ConstantParameter<double> const c{"a", {0., 0., 1.}};
    ParameterLib::CoordinateSystem const coordinate_system{b, c, a};

    std::string m =
        "<medium>"
        "<phases><phase><type>AqueousLiquid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.0</value>"
        "  </property> "
        "</properties>"
        "</phase>"
        "<phase><type>Solid</type>"
        "<properties>"
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>Constant</type>"
        "    <value>0.923 0.531 0.89</value>"
        "  </property> "
        "</properties>"
        "</phase></phases>"
        "<properties>"
        "  <property>"
        "    <name>porosity</name>"
        "    <type>Constant</type>"
        "    <value>0.12</value>"
        "  </property> "
        "  <property>"
        "    <name>thermal_conductivity</name>"
        "    <type>EffectiveThermalConductivityPorosityMixing</type>"
        "  </property> "
        "</properties>"
        "</medium>";

    auto const& medium =
        Tests::createTestMaterial(std::move(m), 3, &coordinate_system);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = 1.0;
    variable_array[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
        0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond(0, 0), 0.467280, 1.e-10);
    ASSERT_NEAR(eff_th_cond(1, 1), 0.7832, 1.e-10);
    ASSERT_NEAR(eff_th_cond(2, 2), 0.81224, 1.e-10);
}
