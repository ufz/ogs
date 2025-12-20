// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

TEST(BaseLibStringTools, TryParseVector_ParsesSimpleDoubles)
{
    std::size_t bad = 999;
    auto v = BaseLib::tryParseVector<double>("1 2 3", &bad);
    ASSERT_TRUE(v.has_value());
    EXPECT_THAT(*v, testing::ElementsAre(1.0, 2.0, 3.0));
    EXPECT_EQ(bad, 999u);  // unchanged
}

TEST(BaseLibStringTools, TryParseVector_WhitespaceSignsAndScientific)
{
    std::size_t bad = 0;
    auto v = BaseLib::tryParseVector<double>("  -1\t 2.5\n 3.0e-1  ", &bad);
    ASSERT_TRUE(v.has_value());
    ASSERT_EQ(v->size(), 3u);
    EXPECT_NEAR((*v)[0], -1.0, 1e-15);
    EXPECT_NEAR((*v)[1], 2.5, 1e-15);
    EXPECT_NEAR((*v)[2], 0.3, 1e-15);
}

TEST(BaseLibStringTools, TryParseVector_EmptyStringYieldsEmptyVector)
{
    std::size_t bad = 999;
    auto v = BaseLib::tryParseVector<double>("", &bad);
    ASSERT_TRUE(v.has_value());
    EXPECT_TRUE(v->empty());
    EXPECT_EQ(bad, 999u);  // unchanged
}

TEST(BaseLibStringTools, TryParseVector_OnlyWhitespaceYieldsEmptyVector)
{
    std::size_t bad = 7;
    auto v = BaseLib::tryParseVector<double>(" \t\n ", &bad);
    ASSERT_TRUE(v.has_value());
    EXPECT_TRUE(v->empty());
    EXPECT_EQ(bad, 7u);
}

TEST(BaseLibStringTools, TryParseVector_TrailingGarbageFailsAndReportsIndex)
{
    std::size_t bad = 0;
    auto v = BaseLib::tryParseVector<double>("1 2 3 xyz", &bad);
    EXPECT_FALSE(v.has_value());
    EXPECT_EQ(bad, 4u);  // error at token 4
}

TEST(BaseLibStringTools, TryParseVector_LeadingGarbageFailsAtToken1)
{
    std::size_t bad = 0;
    auto v = BaseLib::tryParseVector<double>("foo 1 2", &bad);
    EXPECT_FALSE(v.has_value());
    EXPECT_EQ(bad, 1u);  // error at token 1
}

TEST(BaseLibStringTools, TryParseVector_MixedNumberThenGarbageFailsAt2)
{
    std::size_t bad = 0;
    auto v = BaseLib::tryParseVector<double>("1 foo 3", &bad);
    EXPECT_FALSE(v.has_value());
    EXPECT_EQ(bad, 2u);  // error at token 2
}

TEST(BaseLibStringTools, TryParseVector_CommaSeparatedFailsAt2)
{
    std::size_t bad = 0;
    auto v = BaseLib::tryParseVector<double>("1,2 3", &bad);
    EXPECT_FALSE(v.has_value());
    EXPECT_EQ(bad, 2u);  // error at token 2
}
