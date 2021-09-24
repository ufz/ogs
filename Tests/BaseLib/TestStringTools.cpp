/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "BaseLib/StringTools.h"

TEST(BaseLibStringTools, SplitString)
{
    using namespace testing;
    using namespace std;

    // empty string
    EXPECT_THAT(BaseLib::splitString("", ';'), ContainerEq(list<string>{}));

    // no delimiter
    EXPECT_THAT(BaseLib::splitString("a", ';'), ContainerEq(list<string>{"a"}));

    // some delimited values
    EXPECT_THAT(BaseLib::splitString("a,b,c", ','),
                ContainerEq(list<string>{"a", "b", "c"}));

    // leading delimiter
    EXPECT_THAT(BaseLib::splitString(",a,b,c", ','),
                ContainerEq(list<string>{"", "a", "b", "c"}));

    // double delimiters
    EXPECT_THAT(BaseLib::splitString("a,b,,c", ','),
                ContainerEq(list<string>{"a", "b", "", "c"}));

    // trailing delimiters are ignored
    EXPECT_THAT(BaseLib::splitString("a,b,c,", ','),
                ContainerEq(list<string>{"a", "b", "c"}));

    // ... and behave like this if they are on their own:
    EXPECT_THAT(BaseLib::splitString(",", ','), ContainerEq(list<string>{""}));
}
