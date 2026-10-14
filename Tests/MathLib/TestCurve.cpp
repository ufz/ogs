// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
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
    return MathLib::createCurve(conf.getConfigSubtree("curve"));
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
