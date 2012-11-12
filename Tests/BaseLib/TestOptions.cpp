/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file TestOptions.cpp
 *
 * Created on 2012-11-03 by Norihiro Watanabe
 */

// ** INCLUDES **
#include "gtest/gtest.h"

#include "Options.h"

TEST(BaseLib, Options)
{
    BaseLib::Options options;
    BaseLib::Options *child = options.addSubGroup("some group name");
    child->addOption("testStr", "some value");
    child->addOptionAsNum("testInt", 123);
    child->addOptionAsNum("testDouble", 0.11);

    ASSERT_TRUE(0==options.getSubGroup("some false group name"));
    ASSERT_TRUE(0!=options.getSubGroup("some group name"));

    const BaseLib::Options* op = options.getSubGroup("some group name");
    ASSERT_EQ("some value", op->getOption("testStr"));
    ASSERT_EQ(123, op->getOptionAsNum<int>("testInt"));
    ASSERT_EQ(0.11, op->getOptionAsNum<double>("testDouble"));
}

