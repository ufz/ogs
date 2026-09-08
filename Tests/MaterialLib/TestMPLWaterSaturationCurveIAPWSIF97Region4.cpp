// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>
#include <spdlog/fmt/fmt.h>

#include <string>
#include <string_view>

#include "MaterialLib/MPL/Properties/WaterSaturationCurveIAPWSIF97Region4.h"
#include "MaterialLib/PhysicalConstant.h"
#include "NumLib/Exceptions.h"

namespace
{
constexpr double critical_pressure =
    MaterialLib::PhysicalConstant::CriticalPoint::PressureWater;

// The message the check aborts the assembly with, or an empty string if it
// accepted the pressure.
std::string rejectionMessage(double const pressure,
                             std::string_view const quantity)
{
    try
    {
        MaterialPropertyLib::IAPWSIF97Region4::checkPressureInRange(pressure,
                                                                    quantity);
    }
    catch (NumLib::AssemblyException const& e)
    {
        return e.what();
    }
    return {};
}
}  // namespace

// Both bounds belong to the correlation's own domain, so the check must accept
// them: the formulation is defined on the closed interval, and rejecting an
// endpoint would abandon steps that sit exactly on the critical pressure.
TEST(MaterialLibIAPWSIF97Region4, InRangePressureIsAccepted)
{
    using namespace MaterialPropertyLib::IAPWSIF97Region4;

    EXPECT_NO_THROW(
        checkPressureInRange(minimum_saturation_pressure, "a test quantity"));
    EXPECT_NO_THROW(checkPressureInRange(1e6, "a test quantity"));
    EXPECT_NO_THROW(checkPressureInRange(critical_pressure, "a test quantity"));
}

// Outside the range the correlations are extrapolated and no longer describe
// water, so the assembly is abandoned rather than carrying on with a number
// that is not a property of anything, and the time stepper repeats the step.
// The message has to say how far outside the pressure was and which quantity
// was being evaluated, because a retreat that keeps happening is diagnosed
// from these messages alone.
TEST(MaterialLibIAPWSIF97Region4, OutOfRangePressureAbortsTheAssembly)
{
    using namespace MaterialPropertyLib::IAPWSIF97Region4;

    EXPECT_THROW(checkPressureInRange(0.5 * minimum_saturation_pressure,
                                      "the liquid density"),
                 NumLib::AssemblyException);
    EXPECT_THROW(
        checkPressureInRange(2 * critical_pressure, "the vapour enthalpy"),
        NumLib::AssemblyException);

    {
        double const pressure = 0.5 * minimum_saturation_pressure;
        std::string const message =
            rejectionMessage(pressure, "the liquid density");

        EXPECT_NE(std::string::npos,
                  message.find(fmt::format("{:g}", pressure)))
            << "message: " << message;
        EXPECT_NE(std::string::npos, message.find("the liquid density"))
            << "message: " << message;
        EXPECT_NE(
            std::string::npos,
            message.find(fmt::format("{:g}", minimum_saturation_pressure)))
            << "message: " << message;
    }

    {
        double const pressure = 2 * critical_pressure;
        std::string const message =
            rejectionMessage(pressure, "the vapour enthalpy");

        EXPECT_NE(std::string::npos,
                  message.find(fmt::format("{:g}", pressure)))
            << "message: " << message;
        EXPECT_NE(std::string::npos, message.find("the vapour enthalpy"))
            << "message: " << message;
        EXPECT_NE(std::string::npos,
                  message.find(fmt::format("{:g}", critical_pressure)))
            << "message: " << message;
    }
}
