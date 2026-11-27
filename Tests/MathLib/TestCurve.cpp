// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "MathLib/Curve/CreateCurve.h"
#include "MathLib/Curve/Curve.h"
#include "Tests/TestTools.h"

namespace
{
std::unique_ptr<MathLib::Curve> createCurveFromXml(char const xml[])
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    return MathLib::createCurve("test", conf.getConfigSubtree("curve"));
}

void expectSameInterpolant(MathLib::Curve const& expected,
                           MathLib::Curve const& actual,
                           std::vector<double> const& x)
{
    // Sample beyond both ends of the support to cover extrapolation too.
    std::vector<double> query_points = x;
    for (int i = -10; i <= 60; ++i)
    {
        query_points.push_back(i * 0.13);
    }
    for (double const q : query_points)
    {
        EXPECT_EQ(expected.getValue(q), actual.getValue(q)) << q;
        EXPECT_EQ(expected.getDerivative(q), actual.getDerivative(q)) << q;
    }
    EXPECT_EQ(x.front(), actual.getSupportMin());
    EXPECT_EQ(x.back(), actual.getSupportMax());
}
}  // namespace

TEST(MathLibCurve, LinearCreatedFromConfig)
{
    auto const curve = createCurveFromXml(
        "<curve>"
        "   <coords> 0 1 2.5 3</coords>"
        "   <values> 0 1 1.2 4</values>"
        "</curve>");

    EXPECT_NEAR(0.5, curve->getValue(0.5), 1e-14);
    EXPECT_NEAR(1.1, curve->getValue(1.75), 1e-14);
    EXPECT_NEAR(1.0, curve->getDerivative(0.5), 1e-14);
    EXPECT_EQ(0.0, curve->getSupportMin());
    EXPECT_EQ(3.0, curve->getSupportMax());
}

TEST(MathLibCurve, SinglePointCurveIsConstant)
{
    auto const curve = createCurveFromXml(
        "<curve>"
        "   <coords> 5</coords>"
        "   <values> 1</values>"
        "</curve>");

    EXPECT_EQ(1.0, curve->getValue(0.0));
    EXPECT_EQ(1.0, curve->getValue(5.0));
    EXPECT_EQ(1.0, curve->getValue(10.0));

    // In particular at the support point itself, which must not read past the
    // end of the single-element data.
    EXPECT_EQ(0.0, curve->getDerivative(0.0));
    EXPECT_EQ(0.0, curve->getDerivative(5.0));
    EXPECT_EQ(0.0, curve->getDerivative(10.0));

    EXPECT_EQ(5.0, curve->getSupportMin());
    EXPECT_EQ(5.0, curve->getSupportMax());
}

TEST(MathLibCurve, EmptyOrMismatchedDataIsRejected)
{
    EXPECT_THROW(MathLib::Curve(MathLib::PiecewiseLinear{},
                                std::vector<double>{}, std::vector<double>{}),
                 std::runtime_error);
    EXPECT_THROW(
        MathLib::Curve(MathLib::PiecewiseLinear{}, std::vector<double>{0, 1},
                       std::vector<double>{}),
        std::runtime_error);
    EXPECT_THROW(
        MathLib::Curve(MathLib::PiecewiseLinear{}, std::vector<double>{0, 1},
                       std::vector<double>{5}),
        std::runtime_error);
}

TEST(MathLibCurve, CopyKeepsInterpolant)
{
    std::vector<double> const x = {0, 1, 2.5, 3, 5, 6.5};
    std::vector<double> const y = {0, 1, 1.2, 4, 4.1, 2.0};
    MathLib::Curve const curve{MathLib::PiecewiseLinear{}, x, y};
    auto const copy = curve;
    expectSameInterpolant(curve, copy, x);
}

TEST(MathLibCurve, CopyAssignmentKeepsInterpolant)
{
    std::vector<double> const x = {0, 1, 2.5, 3, 5, 6.5};
    std::vector<double> const y = {0, 1, 1.2, 4, 4.1, 2.0};
    MathLib::Curve const source{MathLib::PiecewiseLinear{}, x, y};
    // Different supporting points, so an assignment that kept the target's
    // interpolation or its supporting-point range shows.
    MathLib::Curve target{MathLib::PiecewiseLinear{}, std::vector<double>{1, 2},
                          std::vector<double>{7, 8}};

    target = source;

    expectSameInterpolant(source, target, x);
}

TEST(MathLibCurve, CopyOfMovedFromCurveIsDiagnosed)
{
    MathLib::Curve curve{MathLib::PiecewiseLinear{}, {0, 1, 2}, {0, 1, 4}};
    MathLib::Curve const moved_to{std::move(curve)};
    EXPECT_EQ(4.0, moved_to.getValue(2.0));

    EXPECT_THROW(MathLib::Curve{curve}, std::runtime_error);

    MathLib::Curve target{MathLib::PiecewiseLinear{}, {0, 1}, {0, 1}};
    EXPECT_THROW(target = curve, std::runtime_error);
    // The failed assignment leaves the target unchanged.
    EXPECT_EQ(1.0, target.getValue(1.0));
}

TEST(MathLibCurve, PchipCreatedFromConfig)
{
    auto const curve = createCurveFromXml(
        "<curve>"
        "   <interpolation>pchip</interpolation>"
        "   <coords> 0 1 2.5 3 5 6.5</coords>"
        "   <values> 0 1 1.2 4 4.1 2.0</values>"
        "</curve>");

    // The cubic interpolant, not the piecewise linear one (0.5 at 0.5).
    EXPECT_NEAR(0.59400826446281, curve->getValue(0.5), 1e-12);
    EXPECT_NEAR(1.18801652892562, curve->getDerivative(0.5), 1e-12);
}

TEST(MathLibCurve, InvalidConfigIsRejected)
{
    // End slopes are not allowed for linear interpolation.
    EXPECT_THROW(createCurveFromXml("<curve>"
                                    "   <left_end_slope>0</left_end_slope>"
                                    "   <coords> 0 1 2.5 3</coords>"
                                    "   <values> 0 1 1.2 4</values>"
                                    "</curve>"),
                 std::runtime_error);
    EXPECT_THROW(createCurveFromXml("<curve>"
                                    "   <interpolation>cubic</interpolation>"
                                    "   <coords> 0 1 2.5 3</coords>"
                                    "   <values> 0 1 1.2 4</values>"
                                    "</curve>"),
                 std::runtime_error);
}

TEST(MathLibCurve, LinearInputErrorNamesTheCurve)
{
    try
    {
        createCurveFromXml(
            "<curve>"
            "   <coords> 0 1 1</coords>"
            "   <values> 0 1 2</values>"
            "</curve>");
        FAIL() << "Equal supporting points were accepted.";
    }
    catch (std::runtime_error const& e)
    {
        EXPECT_NE(std::string_view{e.what()}.find("curve 'test'"),
                  std::string_view::npos)
            << e.what();
    }
}

TEST(MathLibCurve, CopyKeepsCubicInterpolant)
{
    std::vector<double> const x = {0, 1, 2.5, 3, 5, 6.5};
    std::vector<double> const y = {0, 1, 1.2, 4, 4.1, 2.0};
    MathLib::Curve const curve{MathLib::Pchip{}, "test", x, y};

    EXPECT_NEAR(0.59400826446281, curve.getValue(0.5), 1e-12);
    auto const copy = curve;
    expectSameInterpolant(curve, copy, x);
}

TEST(MathLibCurve, LinearIsDefault)
{
    for (auto const* const xml : {"<curve>"
                                  "   <coords> 0 1 2.5 3</coords>"
                                  "   <values> 0 1 1.2 4</values>"
                                  "</curve>",
                                  "<curve>"
                                  "   <interpolation>linear</interpolation>"
                                  "   <coords> 0 1 2.5 3</coords>"
                                  "   <values> 0 1 1.2 4</values>"
                                  "</curve>"})
    {
        auto const curve = createCurveFromXml(xml);
        EXPECT_NEAR(0.5, curve->getValue(0.5), 1e-14);
        EXPECT_NEAR(1.1, curve->getValue(1.75), 1e-14);
    }
}

TEST(MathLibCurve, PchipEndSlopesFromConfig)
{
    auto const curve = createCurveFromXml(
        "<curve>"
        "   <interpolation>pchip</interpolation>"
        "   <left_end_slope>0</left_end_slope>"
        "   <right_end_slope>-0.5</right_end_slope>"
        "   <coords> 0 1 2.5 3 5 6.5</coords>"
        "   <values> 0 1 1.2 4 4.1 2.0</values>"
        "</curve>");

    EXPECT_NEAR(0.0, curve->getDerivative(0), 1e-14);
    EXPECT_NEAR(-0.5, curve->getDerivative(6.5), 1e-14);
}
