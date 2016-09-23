/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "ProcessLib/Parameter/PropertyIndexParameter.h"

// make useful line numbers appear in the output of gtest
#define EXPECT_ERR_WARN(cbs, error, warning) do  { \
        if (error)   EXPECT_TRUE((cbs).get_error());   else EXPECT_FALSE((cbs).get_error()); \
        if (warning) EXPECT_TRUE((cbs).get_warning()); else EXPECT_FALSE((cbs).get_warning()); \
        (cbs).reset(); \
    } while(false)

namespace
{

// Exception thrown by the error callback of the class below
class Exc {};

// class that provides callback functions used with ConfigTree
class Callbacks
{
public:
    BaseLib::ConfigTree::Callback
    get_error_cb()
    {
        return [this](std::string const& filename, std::string const& path,
                      std::string const& message)
        {
            (void) path; (void) message;

            // check that filename is passed around properly, especially with
            // move construction/assignment
            EXPECT_EQ("FILENAME", filename);

            DBUG("error <%s> : %s", path.c_str(), message.c_str());
            _error = true;
            throw Exc(); // throw in order to stop normal execution
        };
    }

    BaseLib::ConfigTree::Callback
    get_warning_cb()
    {
        return [this](std::string const& filename, std::string const& path,
                      std::string const& message)
        {
            (void) path; (void) message;

            // check that filename is passed around properly, especially with
            // move construction/assignment
            EXPECT_EQ("FILENAME", filename);

            DBUG("warning <%s> : %s", path.c_str(), message.c_str());
            _warning = true;
        };
    }

    bool get_error()   const { return _error; }
    bool get_warning() const { return _warning; }
    void reset() { _error = false; _warning = false; }

private:
    bool _error = false;
    bool _warning = false;
};

BaseLib::ConfigTree
makeConfigTree(boost::property_tree::ptree const& ptree, Callbacks& cbs)
{
    return BaseLib::ConfigTree(ptree, "FILENAME",
                                  cbs.get_error_cb(), cbs.get_warning_cb());
}

} // no named namespace

TEST(ProcessLib_Parameter, PropertyIndexParameterElement)
{
    const char xml[] =
            "<parameter>"
            "<type>PropertyIndex</type>"
            "<property_name>MaterialIDs</property_name>"
            "<index_values><index>0</index><value>0</value></index_values>"
            "<index_values><index>1</index><value>100</value></index_values>"
            "<index_values><index>3</index><value>300</value></index_values>"
            "</parameter>";
    auto const ptree = readXml(xml);

    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::MeshGenerator::generateLineMesh(4u, 1.0));
    std::vector<int> mat_ids({0,1,2,3});
    MeshLib::addPropertyToMesh(*mesh, "MaterialIDs", MeshLib::MeshItemType::Cell, 1, mat_ids);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(ptree, cbs);
        std::unique_ptr<ProcessLib::ParameterBase> parameter_base = ProcessLib::createPropertyIndexParameter(conf.getConfigSubtree("parameter"), *mesh);
        EXPECT_ERR_WARN(cbs, false, false);

        auto parameter = dynamic_cast<ProcessLib::Parameter<double>*>(parameter_base.get());
        double t = 0;
        ProcessLib::SpatialPosition x;
        x.setElementID(0);
        ASSERT_EQ(0.0, (*parameter)(t, x)[0]);
        x.setElementID(1);
        ASSERT_EQ(100.0, (*parameter)(t, x)[0]);
        x.setElementID(2);
        ASSERT_TRUE((*parameter)(t, x).empty());
        x.setElementID(3);
        ASSERT_EQ(300.0, (*parameter)(t, x)[0]);


    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(ProcessLib_Parameter, PropertyIndexParameterNode)
{
    const char xml[] =
            "<parameter>"
            "<type>PropertyIndex</type>"
            "<property_name>PointGroupIDs</property_name>"
            "<index_values><index>0</index><value>0</value></index_values>"
            "<index_values><index>1</index><value>100</value></index_values>"
            "<index_values><index>3</index><value>300</value></index_values>"
            "</parameter>";
    auto const ptree = readXml(xml);

    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::MeshGenerator::generateLineMesh(4u, 1.0));
    std::vector<int> group_ids({0,1,2,3,4});
    MeshLib::addPropertyToMesh(*mesh, "PointGroupIDs", MeshLib::MeshItemType::Node, 1, group_ids);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(ptree, cbs);
        std::unique_ptr<ProcessLib::ParameterBase> parameter_base = ProcessLib::createPropertyIndexParameter(conf.getConfigSubtree("parameter"), *mesh);
        EXPECT_ERR_WARN(cbs, false, false);

        auto parameter = dynamic_cast<ProcessLib::Parameter<double>*>(parameter_base.get());
        double t = 0;
        ProcessLib::SpatialPosition x;
        x.setNodeID(0);
        ASSERT_EQ(0.0, (*parameter)(t, x)[0]);
        x.setNodeID(1);
        ASSERT_EQ(100.0, (*parameter)(t, x)[0]);
        x.setNodeID(2);
        ASSERT_TRUE((*parameter)(t, x).empty());
        x.setNodeID(3);
        ASSERT_EQ(300.0, (*parameter)(t, x)[0]);
        x.setNodeID(4);
        ASSERT_TRUE((*parameter)(t, x).empty());


    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

