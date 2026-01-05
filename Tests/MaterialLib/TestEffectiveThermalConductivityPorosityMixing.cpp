// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/EffectiveThermalConductivityPorosityMixing.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

std::string phase(std::string name, double thermal_conductivity)
{
    std::stringstream m;
    m << "<phase>\n";
    m << "<type>" << name << "</type>";
    m << "<properties>";
    m << "  <property>";
    m << "    <name>thermal_conductivity</name>";
    m << "    <type>Constant</type>";
    m << "    <value>" << thermal_conductivity << "</value>";
    m << "  </property> ";
    m << "</properties>";
    m << "</phase>";
    return m.str();
};

std::string porosity(double phi)
{
    std::stringstream m;
    m << "  <property>";
    m << "    <name>porosity</name>";
    m << "    <type>Constant</type>";
    m << "    <value>" << phi << "</value>";
    m << "  </property> ";
    return m.str();
};

std::string thermal_conductivity()
{
    std::stringstream m;
    m << "  <property> ";
    m << "    <name>thermal_conductivity</name>";
    m << "    <type>EffectiveThermalConductivityPorosityMixing</type>";
    m << "  </property> ";
    return m.str();
}

const std::string material_liquid_solid =
    "<medium><phases>" + phase("AqueousLiquid", 0.123) + phase("Solid", 0.923) +
    "</phases><properties>" + porosity(0.12) + thermal_conductivity() +
    "</properties></medium>";

const std::string material_gas_solid =
    "<medium><phases>" + phase("Gas", 0.123) + phase("Solid", 0.923) +
    "</phases><properties>" + porosity(0.12) + thermal_conductivity() +
    "</properties></medium>";
const std::string material_gas_liquid_solid =
    "<medium><phases>" + phase("Gas", 0.123) + phase("AqueousLiquid", 0.456) +
    phase("Solid", 0.789) + "</phases><properties>" + porosity(0.12) +
    thermal_conductivity() + "</properties></medium>";

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingLiquidSolid1D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_liquid_solid), 1);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 1.0;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<1>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 1, 1.e-10);
}
TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingLiquidSolid2D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_liquid_solid), 2);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 1.0;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<2>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 2, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<Eigen::Matrix<double, 2, 2>>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt);
    EXPECT_PRED_FORMAT2(Tests::EigenIsNear{},
                        decltype(d_eff_th_cond_dT)::Zero(), d_eff_th_cond_dT);
}
TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingLiquidSolid3D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_liquid_solid), 3);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 1.0;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 3, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<Eigen::Matrix<double, 3, 3>>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt);
    EXPECT_PRED_FORMAT2(Tests::EigenIsNear{},
                        decltype(d_eff_th_cond_dT)::Zero(), d_eff_th_cond_dT);
}

TEST(MaterialPropertyLib, EffectiveThermalConductivityPorosityMixingGasSolid1D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_gas_solid), 1);
    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 0.0;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<1>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 1, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<double>(variable_array,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, time, dt);
    ASSERT_NEAR(d_eff_th_cond_dT, 0.0, 1.e-10);
}

TEST(MaterialPropertyLib, EffectiveThermalConductivityPorosityMixingGasSolid2D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_gas_solid), 2);
    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 0.0;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<2>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 2, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<Eigen::Matrix<double, 2, 2>>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt);
    EXPECT_PRED_FORMAT2(Tests::EigenIsNear{},
                        decltype(d_eff_th_cond_dT)::Zero(), d_eff_th_cond_dT);
}

TEST(MaterialPropertyLib, EffectiveThermalConductivityPorosityMixingGasSolid3D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_gas_solid), 3);
    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 0.0;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.827 * 3, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<Eigen::Matrix<double, 3, 3>>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt);
    EXPECT_PRED_FORMAT2(Tests::EigenIsNear{},
                        decltype(d_eff_th_cond_dT)::Zero(), d_eff_th_cond_dT);
}

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingGasLiquidSolid1D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_gas_liquid_solid), 1);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 0.3;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<1>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.721068 * 1, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<double>(variable_array,
                                     MaterialPropertyLib::Variable::temperature,
                                     pos, time, dt);
    ASSERT_NEAR(d_eff_th_cond_dT, 0.0, 1.e-10);
}

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingGasLiquidSolid2D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_gas_liquid_solid), 2);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 0.3;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<2>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.721068 * 2, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<Eigen::Matrix<double, 2, 2>>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt);
    EXPECT_PRED_FORMAT2(Tests::EigenIsNear{},
                        decltype(d_eff_th_cond_dT)::Zero(), d_eff_th_cond_dT);
}

TEST(MaterialPropertyLib,
     EffectiveThermalConductivityPorosityMixingGasLiquidSolid3D)
{
    auto const& medium =
        Tests::createTestMaterial(std::move(material_gas_liquid_solid), 3);

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.liquid_saturation = 0.3;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond.trace(), 0.721068 * 3, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<Eigen::Matrix<double, 3, 3>>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt);
    EXPECT_PRED_FORMAT2(Tests::EigenIsNear{},
                        decltype(d_eff_th_cond_dT)::Zero(), d_eff_th_cond_dT);
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

    variable_array.liquid_saturation = 1.0;
    variable_array.porosity = 0.12;
    auto const eff_th_cond = MaterialPropertyLib::formEigenTensor<3>(
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, time, dt));
    ASSERT_NEAR(eff_th_cond(0, 0), 0.7832, 1.e-10);
    ASSERT_NEAR(eff_th_cond(1, 1), 0.81224, 1.e-10);
    ASSERT_NEAR(eff_th_cond(2, 2), 0.467280, 1.e-10);

    // Test derivative w.r.t. temperature - should be zero for constant
    // properties
    auto const d_eff_th_cond_dT =
        medium
            ->property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<Eigen::Matrix<double, 3, 3>>(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt);
    EXPECT_PRED_FORMAT2(Tests::EigenIsNear{},
                        decltype(d_eff_th_cond_dT)::Zero(), d_eff_th_cond_dT);
}
