// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <optional>
#include <stdexcept>

#include "ProcessLib/Utils/CouplingScheme.h"

TEST(ProcessLib, CouplingSchemeUnsetIsMonolithic)
{
    EXPECT_TRUE(ProcessLib::useMonolithicScheme(std::nullopt, "Test"));
}

TEST(ProcessLib, CouplingSchemeMonolithic)
{
    EXPECT_TRUE(
        ProcessLib::useMonolithicScheme(std::string{"monolithic"}, "Test"));
}

TEST(ProcessLib, CouplingSchemeDeprecatedMonolithicScheme)
{
    EXPECT_TRUE(ProcessLib::useMonolithicScheme(
        std::string{"monolithic_scheme"}, "Test"));
}

TEST(ProcessLib, CouplingSchemeStaggered)
{
    EXPECT_FALSE(
        ProcessLib::useMonolithicScheme(std::string{"staggered"}, "Test"));
}

TEST(ProcessLib, CouplingSchemeUnknownValueIsFatal)
{
    EXPECT_THROW(
        ProcessLib::useMonolithicScheme(std::string{"Staggered"}, "Test"),
        std::runtime_error);
}
