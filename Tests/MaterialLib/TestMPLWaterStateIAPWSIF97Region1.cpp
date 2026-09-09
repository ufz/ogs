// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>
#include <spdlog/fmt/fmt.h>

#include <string>
#include <string_view>

#include "MaterialLib/MPL/Properties/WaterStateIAPWSIF97Region1.h"
#include "NumLib/Exceptions.h"

namespace
{
// The message the check aborts the assembly with, or an empty string if it
// accepted the state.
std::string rejectionMessage(double const pressure, double const temperature,
                             std::string_view const quantity)
{
    try
    {
        MaterialPropertyLib::IAPWSIF97Region1::checkStateInRange(
            pressure, temperature, quantity);
    }
    catch (NumLib::AssemblyException const& e)
    {
        return e.what();
    }
    return {};
}
}  // namespace

// Both bounds belong to region 1 itself, so a state sitting exactly on one of
// them is compressed liquid and must not be abandoned. The pressures above the
// critical pressure are the ones the wellbore process reaches this check with
// at all, since below it the saturation line answers instead.
TEST(MaterialLibIAPWSIF97Region1, InRangeStateIsAccepted)
{
    using namespace MaterialPropertyLib::IAPWSIF97Region1;

    EXPECT_NO_THROW(checkStateInRange(30e6, 500., "a test quantity"));
    EXPECT_NO_THROW(
        checkStateInRange(maximum_pressure, 500., "a test quantity"));
    EXPECT_NO_THROW(
        checkStateInRange(30e6, maximum_temperature, "a test quantity"));
    EXPECT_NO_THROW(checkStateInRange(maximum_pressure, maximum_temperature,
                                      "a test quantity"));
}

// Beyond either bound the correlations describe region 2 or region 3 rather
// than the compressed liquid, and they do not say so: they return plausible
// looking values instead. The caller has already established that no
// saturation state exists at this pressure, so there is nothing to fall back
// on and the assembly is abandoned, which lets the time stepper repeat the
// step with a smaller step size.
TEST(MaterialLibIAPWSIF97Region1, OutOfRangeStateAbortsTheAssembly)
{
    using namespace MaterialPropertyLib::IAPWSIF97Region1;

    EXPECT_THROW(
        checkStateInRange(2 * maximum_pressure, 500., "a test quantity"),
        NumLib::AssemblyException);
    EXPECT_THROW(
        checkStateInRange(30e6, 2 * maximum_temperature, "a test quantity"),
        NumLib::AssemblyException);
}

// A retreat that keeps happening is diagnosed from these messages alone, so
// the message carries both parts of the state, both bounds, and what was being
// evaluated.
TEST(MaterialLibIAPWSIF97Region1, MessageReportsTheStateAndTheBounds)
{
    using namespace MaterialPropertyLib::IAPWSIF97Region1;

    std::string const message =
        rejectionMessage(2 * maximum_pressure, 700., "the drift-flux closure");

    EXPECT_NE(std::string::npos, message.find(fmt::format("{:g}", 700.)))
        << "message: " << message;
    EXPECT_NE(std::string::npos,
              message.find(fmt::format("{:g}", 2 * maximum_pressure)))
        << "message: " << message;
    EXPECT_NE(std::string::npos,
              message.find(fmt::format("{:g}", maximum_pressure)))
        << "message: " << message;
    EXPECT_NE(std::string::npos,
              message.find(fmt::format("{:g}", maximum_temperature)))
        << "message: " << message;
    EXPECT_NE(std::string::npos, message.find("the drift-flux closure"))
        << "message: " << message;
}
