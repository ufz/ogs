/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include <boost/property_tree/xml_parser.hpp>
#include <numeric>
#include <sstream>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "Tests/TestTools.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "ProcessLib/Parameter/GroupBasedParameter.h"

using namespace ProcessLib;

std::unique_ptr<Parameter<double>> constructParameterFromString(
    std::string const& xml, std::vector<MeshLib::Mesh*> const meshes,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves = {})
{
    auto const xml_ptree = readXml(xml.c_str());
    BaseLib::ConfigTree config_tree(xml_ptree, "", BaseLib::ConfigTree::onerror,
                                    BaseLib::ConfigTree::onwarning);
    auto parameter_base = createParameter(config_tree, meshes, curves);
    return std::unique_ptr<Parameter<double>>(
        static_cast<Parameter<double>*>(parameter_base.release()));
}

struct ProcessLibParameter : public ::testing::Test
{
    void SetUp() override
    {
        // A mesh with four elements, five points.
        mesh.reset(MeshLib::MeshGenerator::generateLineMesh(4u, 1.0));
    }
    std::unique_ptr<MeshLib::Mesh> mesh;
};

TEST_F(ProcessLibParameter, GroupBasedParameterElement)
{
    std::vector<int> mat_ids({0, 1, 2, 3});
    MeshLib::addPropertyToMesh(*mesh, "MaterialIDs",
                               MeshLib::MeshItemType::Cell, 1, mat_ids);

    auto parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>Group</type>"
        "<group_id_property>MaterialIDs</group_id_property>"
        "<index_values><index>0</index><value>0</value></index_values>"
        "<index_values><index>1</index><value>100</value></index_values>"
        "<index_values><index>3</index><value>300</value></index_values>",
        {mesh.get()});

    double t = 0;
    ProcessLib::SpatialPosition x;
    x.setElementID(0);
    ASSERT_EQ(0.0, (*parameter)(t, x)[0]);
    x.setElementID(1);
    ASSERT_EQ(100.0, (*parameter)(t, x)[0]);
    x.setElementID(2);
    ASSERT_ANY_THROW((*parameter)(t, x));
    x.setElementID(3);
    ASSERT_EQ(300.0, (*parameter)(t, x)[0]);

}

TEST_F(ProcessLibParameter, GroupBasedParameterNode)
{
    std::vector<int> group_ids({0, 1, 2, 3, 4});
    MeshLib::addPropertyToMesh(*mesh, "PointGroupIDs",
                               MeshLib::MeshItemType::Node, 1, group_ids);

    auto parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>Group</type>"
        "<group_id_property>PointGroupIDs</group_id_property>"
        "<index_values><index>0</index><value>0</value></index_values>"
        "<index_values><index>1</index><value>100</value></index_values>"
        "<index_values><index>3</index><value>300</value></index_values>",
        {mesh.get()});

    double t = 0;
    ProcessLib::SpatialPosition x;
    x.setNodeID(0);
    ASSERT_EQ(0.0, (*parameter)(t, x)[0]);
    x.setNodeID(1);
    ASSERT_EQ(100.0, (*parameter)(t, x)[0]);
    x.setNodeID(2);
    ASSERT_ANY_THROW((*parameter)(t, x));
    x.setNodeID(3);
    ASSERT_EQ(300.0, (*parameter)(t, x)[0]);
    x.setNodeID(4);
    ASSERT_ANY_THROW((*parameter)(t, x));
}

