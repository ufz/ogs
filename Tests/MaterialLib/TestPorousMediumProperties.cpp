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

namespace
{
using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

//----------------------------------------------------------------------------
// Test density models.
Storage* createTestStorageModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("storage");
    return MaterialLib::PorousMedium::createStorageModel(&sub_config);
}

TEST(Material, checkConstantStrage)
{
    const char xml[] =
        "<storage>"
        "   <type>constant</type>"
        "   <value> 1.e-4 </value> "
        "</storage>";
    std::unique_ptr<Storage> eta(createTestStorageModel(xml));

    ASSERT_EQ(1.e-4, eta->getValue());
}

}
