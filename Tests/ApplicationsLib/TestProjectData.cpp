/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>

#include "Applications/ApplicationsLib/ProjectData.h"

TEST(ApplicationsLib_ProjectData_SplitIntegerList, EmptyList)
{
    ASSERT_TRUE(splitMaterialIdString("").empty());
}

TEST(ApplicationsLib_ProjectData_SplitIntegerList, SingleInt)
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

TEST(ApplicationsLib_ProjectData_SplitIntegerList, SingleIntFail)
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

    // range exceeeded
    EXPECT_THROW(
        splitMaterialIdString("1234567890123456789012345678901234567890"),
        std::runtime_error);
}

TEST(ApplicationsLib_ProjectData_SplitIntegerList, IntList)
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

TEST(ApplicationsLib_ProjectData_SplitIntegerList, IntListFail)
{
    // empty element
    EXPECT_THROW(splitMaterialIdString("5,,6"), std::runtime_error);

    // leading comma
    EXPECT_THROW(splitMaterialIdString(",40"), std::runtime_error);

    // missing comma
    EXPECT_THROW(splitMaterialIdString("12   20"), std::runtime_error);

    // wrong number in the list
    EXPECT_THROW(splitMaterialIdString("1,2,x,5"), std::runtime_error);
}
