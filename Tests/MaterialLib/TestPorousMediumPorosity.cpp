/*!
   \file  TestPorousMediumPorosity.cpp
   \brief Test the classes for porosity models.

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include <gtest/gtest.h>

#include <memory>

#include "Tests/TestTools.h"

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/PorousMedium/Porosity/ConstantPorosity.h"
#include "MaterialLib/PorousMedium/Porosity/createPorosityModel.h"

using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

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
        "   <type>Constant</type>"
        "   <value> 0.2 </value> "
        "</porosity>";
    auto const n = createTestPorosityModel(xml);

    const double variable = 0.;
    const double temperature = 0.;
    ASSERT_EQ(0.2, n->getValue(variable, temperature));
}
