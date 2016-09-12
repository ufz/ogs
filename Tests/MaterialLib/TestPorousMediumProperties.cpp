/*!
   \file  TestPorousMediumProperties.cpp
   \brief Test classes for fluid properties.

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include <gtest/gtest.h>

#include <memory>

#include "TestTools.h"

#include "MaterialLib/PorousMedium/Storage/createStorageModel.h"
#include "MaterialLib/PorousMedium/Porosity/createPorosityModel.h"

namespace
{
using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

//----------------------------------------------------------------------------
// Test storage models.
std::unique_ptr<Storage> createTestStorageModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("storage");
    return MaterialLib::PorousMedium::createStorageModel(sub_config);
}

TEST(Material, checkConstantStorage)
{
    const char xml[] =
        "<storage>"
        "   <type>constant</type>"
        "   <value> 1.e-4 </value> "
        "</storage>";
    auto const eta = createTestStorageModel(xml);

    ASSERT_EQ(1.e-4, eta->getValue());
}

//----------------------------------------------------------------------------
// Test porosity models.
std::unique_ptr<Porosity> createTestPorosityModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("porosity");
    return MaterialLib::PorousMedium::createPorosityModel(sub_config);
}

TEST(Material, checkConstantPorosity)
{
    const char xml[] =
        "<porosity>"
        "   <type>constant</type>"
        "   <value> 0.2 </value> "
        "</porosity>";
    auto const n = createTestPorosityModel(xml);

    ASSERT_EQ(0.2, n->getValue());
}
}
