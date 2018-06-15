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
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"
#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Parameter/GroupBasedParameter.h"
#include "ProcessLib/Parameter/CurveScaledParameter.h"

using namespace ProcessLib;

std::unique_ptr<Parameter<double>> constructParameterFromString(
    std::string const& xml,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
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
        meshes.emplace_back(MeshLib::MeshGenerator::generateLineMesh(4u, 1.0));
    }
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
};

TEST_F(ProcessLibParameter, GroupBasedParameterElement)
{
    std::vector<int> mat_ids({0, 1, 2, 3});
    MeshLib::addPropertyToMesh(*meshes[0], "MaterialIDs",
                               MeshLib::MeshItemType::Cell, 1, mat_ids);

    auto parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>Group</type>"
        "<group_id_property>MaterialIDs</group_id_property>"
        "<index_values><index>0</index><value>0</value></index_values>"
        "<index_values><index>1</index><value>100</value></index_values>"
        "<index_values><index>3</index><value>300</value></index_values>",
        meshes);

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
    MeshLib::addPropertyToMesh(*meshes[0], "PointGroupIDs",
                               MeshLib::MeshItemType::Node, 1, group_ids);

    auto parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>Group</type>"
        "<group_id_property>PointGroupIDs</group_id_property>"
        "<index_values><index>0</index><value>0</value></index_values>"
        "<index_values><index>1</index><value>100</value></index_values>"
        "<index_values><index>3</index><value>300</value></index_values>",
        meshes);

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

bool testNodalValuesOfElement(
    std::vector<MeshLib::Element*> const& elements,
    std::function<double(MeshLib::Element* const e,
                         std::size_t const local_node_id)>
        get_expected_value,
    Parameter<double> const& parameter,
    double const t)
{
    return std::all_of(
        begin(elements), end(elements), [&](MeshLib::Element* const e) {
            // scalar values, 2 nodes.
            Eigen::Matrix<double, 2, 1> const nodal_values =
                parameter.getNodalValuesOnElement(*e, t);
            if (nodal_values.cols() != 1)
            {
                ERR("Expected scalar nodal values.");
                return false;
            }
            if (nodal_values.rows() != static_cast<int>(e->getNumberOfNodes()))
            {
                ERR("Expected equal number of element nodes and the nodal "
                    "values.");
                return false;
            }
            for (std::size_t i = 0; i < e->getNumberOfNodes(); ++i)
            {
                double const expected_value = get_expected_value(e, i);
                if (expected_value != nodal_values(i, 0))
                {
                    ERR("Mismatch for element %d, node %d; Expected %g, got "
                        "%g.",
                        e->getID(), i, expected_value, nodal_values(i, 0));
                    return false;
                }
            }
            return true;
        });
}

// For all elements all nodes have a constant value.
TEST_F(ProcessLibParameter, GetNodalValuesOnElement_constant)
{
    auto const parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>Constant</type>"
        "<value>42.23</value>",
        meshes);

    double const t = 0;
    auto expected_value = [](MeshLib::Element* const /*e*/,
                              std::size_t const /*local_node_id*/) {
        return 42.23;
    };

    ASSERT_TRUE(testNodalValuesOfElement(meshes[0]->getElements(),
                                         expected_value, *parameter, t));
}

TEST_F(ProcessLibParameter, GetNodalValuesOnElement_node)
{
    std::vector<double> node_ids({0, 1, 2, 3, 4});
    MeshLib::addPropertyToMesh(*meshes[0], "NodeIDs",
                               MeshLib::MeshItemType::Node, 1, node_ids);

    auto const parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>MeshNode</type>"
        "<field_name>NodeIDs</field_name>",
        meshes);

    double const t = 0;

    // For all elements all nodes have the value of the node id.
    auto expected_value = [](MeshLib::Element* const e,
                              std::size_t const local_node_id) {
        return static_cast<double>(e->getNode(local_node_id)->getID());
    };

    ASSERT_TRUE(testNodalValuesOfElement(meshes[0]->getElements(),
                                         expected_value, *parameter, t));
}

TEST_F(ProcessLibParameter, GetNodalValuesOnElement_element)
{
    std::vector<double> element_ids({0, 1, 2, 3});
    MeshLib::addPropertyToMesh(*meshes[0], "ElementIDs",
                               MeshLib::MeshItemType::Cell, 1, element_ids);

    auto const parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>MeshElement</type>"
        "<field_name>ElementIDs</field_name>",
        meshes);

    double const t = 0;

    // For all elements all nodes have the value of the element id.
    auto expected_value = [](MeshLib::Element* const e,
                             std::size_t const /*local_node_id*/) {
        return static_cast<double>(e->getID());
    };

    ASSERT_TRUE(testNodalValuesOfElement(meshes[0]->getElements(),
                                         expected_value, *parameter, t));
}

TEST_F(ProcessLibParameter, GetNodalValuesOnElement_curveScaledNode)
{
    std::vector<double> node_ids({0, 1, 2, 3, 4});
    MeshLib::addPropertyToMesh(*meshes[0], "NodeIDs",
                               MeshLib::MeshItemType::Node, 1, node_ids);

    std::vector<std::unique_ptr<ParameterBase>> parameters;
    parameters.emplace_back(
        constructParameterFromString("<name>NodeIDs</name>"
                                     "<type>MeshNode</type>"
                                     "<field_name>NodeIDs</field_name>",
                                     meshes));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 1}, std::vector<double>{0, 1}, true);

    auto const parameter = constructParameterFromString(
        "<name>parameter</name>"
        "<type>CurveScaled</type>"
        "<curve>linear_curve</curve>"
        "<parameter>NodeIDs</parameter>",
        meshes, curves);

    parameter->initialize(parameters);

    double const t = 0.5;

    // For all elements all nodes have the value of the node id times the time.
    auto expected_value = [&t](MeshLib::Element* const e,
                              std::size_t const local_node_id) {
        return static_cast<double>(e->getNode(local_node_id)->getID()) * t;
    };

    ASSERT_TRUE(testNodalValuesOfElement(meshes[0]->getElements(),
                                         expected_value, *parameter, t));
}
