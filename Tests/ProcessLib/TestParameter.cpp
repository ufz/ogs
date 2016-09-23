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

#include "ProcessLib/Parameter/MaterialIDParameter.h"


TEST(ProcessLib_Parameter, MaterialIDParameter)
{
    const char xml[] =
            "<parameter>"
            "<type>MaterialID</type>"
            "<name>ParameterX</name>"
            "<property_name>MaterialIDs</property_name>"
            "<index_values><index>0</index><value>0</value></index_values>"
            "<index_values><index>1</index><value>100</value></index_values>"
            "<index_values><index>3</index><value>300</value></index_values>"
            "</parameter>";
    auto const ptree = readXml(xml);

    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::MeshGenerator::generateLineMesh(4u, 1.0));
    std::vector<int> mat_ids({0,1,2,3});
    MeshLib::addPropertyToMesh(*mesh, "MaterialIDs", MeshLib::MeshItemType::Cell, 1, mat_ids);

    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror, BaseLib::ConfigTree::onwarning);
    std::unique_ptr<ProcessLib::ParameterBase> parameter_base
            = ProcessLib::createMaterialIDParameter(conf.getConfigSubtree("parameter"), *mesh);

    auto parameter = dynamic_cast<ProcessLib::Parameter<double>*>(parameter_base.get());
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
