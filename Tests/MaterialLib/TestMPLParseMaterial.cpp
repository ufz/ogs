/**
 * \author Norbert Grunwald
 * \date   Sep 22, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <memory>

#include "Tests/TestTools.h"

#include "MaterialLib/MPL/mpMedium.h"

namespace MPL = MaterialPropertyLib;

//----------------------------------------------------------------------------
// Test density models.
MPL::Medium createTestMaterial(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config =
            conf.getConfigSubtree("media").getConfigSubtree("medium");

    return MPL::Medium(sub_config);
}

TEST(Material, parseMaterials)
{
    const char xml[] =
            "<media>"
            "  <medium>"
            "    <phases>"
            "      <phase>"
            "        <components>"
            "          <component>"
            "            <properties>"
            "              <property>"
            "                <name>density</name>"
            "                <type>constant</type>"
            "                <value>1234.5</value>"
            "              </property>"
            "            </properties>"
            "          </component>"
            "        </components>"
            "      </phase>"
            "    </phases>"
            "  </medium>"
            "</media>";


    auto m = createTestMaterial(xml);

    const auto mediumName = MPL::getString(m.property(MPL::name));


    ASSERT_EQ(mediumName, "no_name");
//    ArrayType dummy;
//    ASSERT_EQ(998.0, rho->getValue(dummy));
//    ASSERT_EQ(
//        0.0,
//        rho->getdValue(dummy, MaterialLib::Fluid::PropertyVariableType::T));
}

