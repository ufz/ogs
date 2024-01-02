/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <map>
#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"
#include "GeoLib/Raster.h"
#include "InfoLib/TestInfo.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/ComponentTransport/CreateLookupTable.h"
#include "ProcessLib/ComponentTransport/LookupTable.h"
#include "ProcessLib/ProcessVariable.h"
#include "Tests/TestTools.h"

using namespace ProcessLib::ComponentTransport;

namespace
{
std::unique_ptr<LookupTable> createTestLookupTable(const char xml[])
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    {
        std::size_t const n_elements = 4;
        auto mesh = std::unique_ptr<MeshLib::Mesh>(
            MeshToolsLib::MeshGenerator::generateRegularQuadMesh(1.0,
                                                                 n_elements));
        meshes.push_back(std::move(mesh));
    }

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;

    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    {
        auto const& para_config = conf.getConfigSubtree("parameter");
        auto p = ParameterLib::createParameter(para_config, meshes, {}, curves);
        parameters.push_back(std::move(p));
    }

    std::vector<ProcessLib::ProcessVariable> variables;
    auto const& process_variables_config =
        conf.getConfigSubtree("process_variables");
    for (auto var_config :
         process_variables_config.getConfigSubtreeList("process_variable"))
    {
        variables.emplace_back(var_config, *meshes[0], meshes, parameters,
                               curves);
    }

    std::vector<
        std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>>
        process_variables;
    {
        for (auto& pv : variables)
        {
            std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>
                per_process_variable;

            per_process_variable.push_back(pv);
            process_variables.push_back(std::move(per_process_variable));
        }
    }

    std::optional<std::string> file_name(TestInfoLib::TestInfo::data_path +
                                         "/FileIO/sparse_table.txt");
    if (!BaseLib::isProjectDirectorySet())
    {
        BaseLib::setProjectDirectory(BaseLib::extractPath(*file_name));
    }

    return createLookupTable(file_name, process_variables);
}

}  // namespace

TEST(ComponentTransport, checkLookupTable)
{
    const char xml[] =
        "<parameter>"
        "     <name>c0</name>"
        "     <type>Constant</type>"
        "     <value>0</value>"
        "</parameter>"
        "<process_variables>"
        "    <process_variable>"
        "        <name>Ni</name>"
        "        <components>1</components>"
        "        <order>1</order>"
        "        <initial_condition>c0</initial_condition>"
        "        <boundary_conditions>"
        "        </boundary_conditions>"
        "    </process_variable>"
        "    <process_variable>"
        "        <name>Np(5)</name>"
        "        <components>1</components>"
        "        <order>1</order>"
        "        <initial_condition>c0</initial_condition>"
        "        <boundary_conditions>"
        "        </boundary_conditions>"
        "    </process_variable>"
        "    <process_variable>"
        "        <name>Th</name>"
        "        <components>1</components>"
        "        <order>1</order>"
        "        <initial_condition>c0</initial_condition>"
        "        <boundary_conditions>"
        "        </boundary_conditions>"
        "    </process_variable>"
        "    <process_variable>"
        "        <name>Ra</name>"
        "        <components>1</components>"
        "        <order>1</order>"
        "        <initial_condition>c0</initial_condition>"
        "        <boundary_conditions>"
        "        </boundary_conditions>"
        "    </process_variable>"
        "</process_variables>";

    auto const lookup_table = createTestLookupTable(xml);
    auto const& tabular_data = lookup_table->tabular_data;

    ASSERT_EQ(1.e-5, tabular_data.at("Ni")[0]);
    ASSERT_EQ(1e-9, tabular_data.at("Ni_prev")[0]);
    ASSERT_EQ(0., tabular_data.at("Ni_new")[0]);

    ASSERT_EQ(1.e-6, tabular_data.at("Np(5)")[0]);
    ASSERT_EQ(1.e-10, tabular_data.at("Np(5)_prev")[0]);
    ASSERT_EQ(0., tabular_data.at("Np(5)_new")[0]);

    ASSERT_EQ(1.e-7, tabular_data.at("Th")[0]);
    ASSERT_EQ(1.e-11, tabular_data.at("Th_prev")[0]);
    ASSERT_EQ(0., tabular_data.at("Th_new")[0]);

    ASSERT_EQ(1.e-8, tabular_data.at("Ra")[0]);
    ASSERT_EQ(1.e-12, tabular_data.at("Ra_prev")[0]);
    ASSERT_EQ(9.977098263915e-13, tabular_data.at("Ra_new")[0]);
}

class ComponentTransport_FieldGetBoundingSeedPointsTests
    : public ::testing::Test
{
protected:
    ProcessLib::ComponentTransport::Field createField(
        std::vector<double> seed_points)
    {
        return ProcessLib::ComponentTransport::Field{{}, seed_points, "", 0};
    }
};

TEST_F(ComponentTransport_FieldGetBoundingSeedPointsTests,
       TestLowerBoundWarning)
{
    auto const field = createField({1.0, 2.0, 3.0, 4.0, 5.0});
    double const value = 0.5;
    auto const pair = field.getBoundingSeedPoints(value);
    EXPECT_DOUBLE_EQ(pair.first, 1.0);
    EXPECT_DOUBLE_EQ(pair.second, 2.0);
}

TEST_F(ComponentTransport_FieldGetBoundingSeedPointsTests,
       TestUpperBoundWarning)
{
    auto const field = createField({1.0, 2.0, 3.0, 4.0, 5.0});
    double const value = 5.5;
    auto pair = field.getBoundingSeedPoints(value);
    EXPECT_DOUBLE_EQ(pair.first, 4.0);
    EXPECT_DOUBLE_EQ(pair.second, 5.0);
}

TEST_F(ComponentTransport_FieldGetBoundingSeedPointsTests, TestMiddleValue)
{
    auto const field = createField({1.0, 2.0, 3.0, 4.0, 5.0});
    double const value = 3.5;
    auto pair = field.getBoundingSeedPoints(value);
    EXPECT_DOUBLE_EQ(pair.first, 3.0);
    EXPECT_DOUBLE_EQ(pair.second, 4.0);
}

TEST_F(ComponentTransport_FieldGetBoundingSeedPointsTests, TestExactValue)
{
    auto const field = createField({1.0, 2.0, 3.0, 4.0, 5.0});
    double const value = 3.0;
    auto pair = field.getBoundingSeedPoints(value);
    EXPECT_DOUBLE_EQ(pair.first, 2.0);
    EXPECT_DOUBLE_EQ(pair.second, 3.0);
}

TEST_F(ComponentTransport_FieldGetBoundingSeedPointsTests, TestTwoElementVector)
{
    auto const field = createField({2.0, 3.0});
    double const value = 3.0;
    auto pair = field.getBoundingSeedPoints(value);
    EXPECT_DOUBLE_EQ(pair.first, 2.0);
    EXPECT_DOUBLE_EQ(pair.second, 3.0);
}

TEST_F(ComponentTransport_FieldGetBoundingSeedPointsTests, TestOneElementVector)
{
    auto const field = createField({3.0});
    double const value = 3.0;
    ASSERT_THROW(field.getBoundingSeedPoints(value), std::runtime_error);
}

TEST_F(ComponentTransport_FieldGetBoundingSeedPointsTests, TestEmptyVector)
{
    auto const field = createField({});
    double const value = 3.0;
    ASSERT_THROW(field.getBoundingSeedPoints(value), std::runtime_error);
}
