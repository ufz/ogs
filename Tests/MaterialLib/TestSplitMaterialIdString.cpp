// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>

#include "MaterialLib/Utils/MediaCreation.h"

using namespace MaterialLib;

TEST(MaterialLib_SplitIntegerList, EmptyList)
{
    ASSERT_TRUE(splitMaterialIdString("").empty());
}

TEST(MaterialLib_SplitIntegerList, SingleInt)
{
    using namespace testing;

    EXPECT_THAT(splitMaterialIdString("5"), ContainerEq(std::vector<int>{5}));

    // leading whitespace
    EXPECT_THAT(splitMaterialIdString("  16"),
                ContainerEq(std::vector<int>{16}));

    // trailing whitespace
    EXPECT_THAT(splitMaterialIdString("23    "),
                ContainerEq(std::vector<int>{23}));

    // negative numbers are OK
    EXPECT_THAT(splitMaterialIdString("-20"),
                ContainerEq(std::vector<int>{-20}));
}

TEST(MaterialLib_SplitIntegerList, SingleIntFail)
{
    // wrong character prefix/suffix
    EXPECT_THROW(splitMaterialIdString("x"), std::runtime_error);
    EXPECT_THROW(splitMaterialIdString(".5"), std::runtime_error);
    EXPECT_THROW(splitMaterialIdString("5?"), std::runtime_error);
    EXPECT_THROW(splitMaterialIdString("7 !"), std::runtime_error);
    EXPECT_THROW(splitMaterialIdString("8   u"), std::runtime_error);

    // hexadecimal numbers are not accepted
    EXPECT_THROW(splitMaterialIdString("0xfa"), std::runtime_error);

    // another integer is not accepted
    EXPECT_THROW(splitMaterialIdString("1 2"), std::runtime_error);

    // range exceeded
    EXPECT_THROW(
        splitMaterialIdString("1234567890123456789012345678901234567890"),
        std::runtime_error);
}

TEST(MaterialLib_SplitIntegerList, IntList)
{
    using namespace testing;

    EXPECT_THAT(splitMaterialIdString("5,6,7"),
                ContainerEq(std::vector<int>{5, 6, 7}));

    // whitespace around comma
    EXPECT_THAT(splitMaterialIdString("9  ,10,  11,12   ,   13"),
                ContainerEq(std::vector<int>{9, 10, 11, 12, 13}));

    // trailing comma is ignored
    EXPECT_THAT(splitMaterialIdString("20, 22, 24,"),
                ContainerEq(std::vector<int>{20, 22, 24}));
}

TEST(MaterialLib_SplitIntegerList, IntListFail)
{
    // only delimiter
    EXPECT_THROW(splitMaterialIdString(","), std::runtime_error);

    // empty element
    EXPECT_THROW(splitMaterialIdString("5,,6"), std::runtime_error);

    // leading comma
    EXPECT_THROW(splitMaterialIdString(",40"), std::runtime_error);

    // missing comma
    EXPECT_THROW(splitMaterialIdString("12   20"), std::runtime_error);

    // wrong number in the list
    EXPECT_THROW(splitMaterialIdString("1,2,x,5"), std::runtime_error);
}

TEST(MaterialLib_SplitIntegerList, RangeSingle)
{
    using namespace testing;

    EXPECT_THAT(splitMaterialIdString("9:13"),
                ContainerEq(std::vector<int>{9, 10, 11, 12, 13}));

    EXPECT_THAT(splitMaterialIdString("-1:2"),
                ContainerEq(std::vector<int>{-1, 0, 1, 2}));

    EXPECT_THAT(splitMaterialIdString("5:5"), ContainerEq(std::vector<int>{5}));

    EXPECT_THAT(splitMaterialIdString("-10:-8"),
                ContainerEq(std::vector<int>{-10, -9, -8}));
}

TEST(MaterialLib_SplitIntegerList, MixedListWithRanges)
{
    using namespace testing;

    EXPECT_THAT(splitMaterialIdString("1:3,5,7:10"),
                ContainerEq(std::vector<int>{1, 2, 3, 5, 7, 8, 9, 10}));

    EXPECT_THAT(splitMaterialIdString("-1:0,5,10:12"),
                ContainerEq(std::vector<int>{-1, 0, 5, 10, 11, 12}));

    EXPECT_THAT(splitMaterialIdString("1:2,4:6,8:10"),
                ContainerEq(std::vector<int>{1, 2, 4, 5, 6, 8, 9, 10}));
}

TEST(MaterialLib_SplitIntegerList, RangeFail)
{
    // end < start
    EXPECT_THROW(splitMaterialIdString("5:3"), std::runtime_error);

    // empty start
    EXPECT_THROW(splitMaterialIdString(":5"), std::runtime_error);

    // empty end
    EXPECT_THROW(splitMaterialIdString("5:"), std::runtime_error);

    // multiple colons
    EXPECT_THROW(splitMaterialIdString("1:2:3"), std::runtime_error);

    // non-integer end
    EXPECT_THROW(splitMaterialIdString("1:x"), std::runtime_error);

    // non-integer start
    EXPECT_THROW(splitMaterialIdString("x:5"), std::runtime_error);

    // whitespace in range
    EXPECT_THROW(splitMaterialIdString("1 :5"), std::runtime_error);

    // invalid character in range
    EXPECT_THROW(splitMaterialIdString("1:?5"), std::runtime_error);
}
