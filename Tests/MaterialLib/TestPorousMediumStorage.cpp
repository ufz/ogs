/*!
   \file
   \brief Test the classes for storage models.

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include <gtest/gtest.h>

#include <memory>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/PorousMedium/Storage/ConstantStorage.h"
#include "MaterialLib/PorousMedium/Storage/createStorageModel.h"
#include "Tests/TestTools.h"

using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

std::unique_ptr<Storage> createTestStorageModel(const char xml[])
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("storage");
    return MaterialLib::PorousMedium::createStorageModel(sub_config);
}

TEST(Material, checkConstantStorage)
{
    const char xml[] =
        "<storage>"
        "   <type>Constant</type>"
        "   <value> 1.e-4 </value> "
        "</storage>";
    auto const eta = createTestStorageModel(xml);

    const double var = 0;
    ASSERT_EQ(1.e-4, eta->getValue(var));
}
