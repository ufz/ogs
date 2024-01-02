/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "ProcessLib/Assembly/MatrixOutput.h"

TEST(ProcessLibLocalMatrixOutput, ParseSetOfSizeTGood)
{
    namespace PAD = ProcessLib::Assembly::detail;

    // empty
    {
        auto const actual = PAD::parseSetOfSizeT("", "");

        EXPECT_TRUE(actual.empty());
    }

    // single value
    {
        auto const actual = PAD::parseSetOfSizeT("40", "");
        std::unordered_set<std::size_t> const expected{40};

        EXPECT_THAT(actual, ::testing::ContainerEq(expected));
    }

    // several values
    {
        auto const actual = PAD::parseSetOfSizeT("2 +3 5 7", "");
        std::unordered_set<std::size_t> const expected{2, 3, 7, 5};

        EXPECT_THAT(actual, ::testing::ContainerEq(expected));
    }
}

TEST(ProcessLibLocalMatrixOutput, ParseSetOfSizeTBad)
{
    namespace PAD = ProcessLib::Assembly::detail;

    // first value bad
    {
        auto const actual = PAD::parseSetOfSizeT("x2 3 5 7", "");

        EXPECT_TRUE(actual.empty());
    }

    // other value bad
    {
        auto const actual = PAD::parseSetOfSizeT("2A 3 5 7", "");
        std::unordered_set<std::size_t> const expected{2};

        EXPECT_THAT(actual, ::testing::ContainerEq(expected));
    }

    // other value bad #2
    {
        auto const actual = PAD::parseSetOfSizeT("2 -3 5 7", "");
        std::unordered_set<std::size_t> const expected{2};

        EXPECT_THAT(actual, ::testing::ContainerEq(expected));
    }

    // other value bad #3
    {
        auto const actual = PAD::parseSetOfSizeT("2 3 5e 7", "");
        std::unordered_set<std::size_t> const expected{2, 3, 5};

        EXPECT_THAT(actual, ::testing::ContainerEq(expected));
    }

    // other value bad #4
    {
        auto const actual = PAD::parseSetOfSizeT("2 3 5 7!", "");
        std::unordered_set<std::size_t> const expected{2, 3, 5, 7};

        EXPECT_THAT(actual, ::testing::ContainerEq(expected));
    }
}

TEST(ProcessLibLocalMatrixOutput, CreateElementPredicate)
{
    namespace PAD = ProcessLib::Assembly::detail;

    // empty
    {
        auto const actual = PAD::createLocalMatrixOutputElementPredicate("");

        EXPECT_FALSE(actual);
    }

    // all true
    {
        auto const actual = PAD::createLocalMatrixOutputElementPredicate("*");

        ASSERT_TRUE(actual);

        for (auto const element_id : {0, 2, 3, 7, 42, 255, 1000001})
        {
            EXPECT_TRUE(actual(element_id));
        }
    }

    // proper element selection
    {
        auto const actual =
            PAD::createLocalMatrixOutputElementPredicate("2 3 5 7");

        ASSERT_TRUE(actual);

        for (auto const element_id : {2, 3, 5, 7})
        {
            EXPECT_TRUE(actual(element_id));
        }
        for (auto const element_id : {0, 1, 4, 10, 10001})
        {
            EXPECT_FALSE(actual(element_id));
        }
    }

    // erroneous element selection
    {
        auto const actual =
            PAD::createLocalMatrixOutputElementPredicate("2 -3 5 7");

        ASSERT_TRUE(actual);

        for (auto const element_id : {2})
        {
            EXPECT_TRUE(actual(element_id));
        }
        for (auto const element_id : {0, 1, 3, 4, 5, 7, 10, 10001})
        {
            EXPECT_FALSE(actual(element_id));
        }
    }

    // completely wrong element selection
    {
        auto const actual =
            PAD::createLocalMatrixOutputElementPredicate("x2 3 5 7");

        ASSERT_FALSE(actual);
    }
}
