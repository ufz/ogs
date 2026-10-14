// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MathLib/InterpolationAlgorithms/SupportingPoints.h"

TEST(MathLibSupportingPoints, SortsValuesAlong)
{
    std::vector<double> x = {3, 1, 2};
    std::vector<double> y = {30, 10, 20};
    MathLib::sortAndCheckSupportingPoints(x, y, 1, "test", false);
    EXPECT_EQ((std::vector<double>{1, 2, 3}), x);
    EXPECT_EQ((std::vector<double>{10, 20, 30}), y);
}

TEST(MathLibSupportingPoints, InvalidInputIsRejected)
{
    auto check = [](std::vector<double> x, std::vector<double> y,
                    std::size_t const min_size, bool const sorted)
    { MathLib::sortAndCheckSupportingPoints(x, y, min_size, "test", sorted); };

    EXPECT_NO_THROW(check({1}, {10}, 1, false));
    // Different numbers of supporting points and values.
    EXPECT_THROW(check({1, 2}, {10}, 1, false), std::runtime_error);
    // Too few supporting points.
    EXPECT_THROW(check({}, {}, 1, false), std::runtime_error);
    EXPECT_THROW(check({1, 2, 3}, {10, 20, 30}, 4, false), std::runtime_error);
    // Equal supporting points, also after sorting.
    EXPECT_THROW(check({1, 2, 1}, {10, 20, 30}, 1, false), std::runtime_error);
    EXPECT_THROW(check({1, 1, 2}, {10, 20, 30}, 1, true), std::runtime_error);
    // Wrongly claimed to be sorted.
    EXPECT_THROW(check({2, 1, 3}, {20, 10, 30}, 1, true), std::runtime_error);
    EXPECT_NO_THROW(check({1, 2, 3}, {10, 20, 30}, 1, true));
    // Non-finite supporting points; NaN would otherwise pass unsorted.
    double const nan = std::numeric_limits<double>::quiet_NaN();
    double const inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(check({3, 1, nan, 2, 0}, {0, 0, 0, 0, 0}, 1, false),
                 std::runtime_error);
    EXPECT_THROW(check({0, 1, nan}, {0, 0, 0}, 1, true), std::runtime_error);
    EXPECT_THROW(check({0, 1, inf}, {0, 0, 0}, 1, false), std::runtime_error);
}

TEST(MathLibSupportingPoints, PiecewiseLinearInterpolationChecksInput)
{
    EXPECT_NO_THROW(MathLib::PiecewiseLinearInterpolation(
        std::vector<double>{1}, std::vector<double>{10}));
    EXPECT_THROW(MathLib::PiecewiseLinearInterpolation(std::vector<double>{},
                                                       std::vector<double>{}),
                 std::runtime_error);
    EXPECT_THROW(MathLib::PiecewiseLinearInterpolation(
                     std::vector<double>{1, 2}, std::vector<double>{10}),
                 std::runtime_error);
    EXPECT_THROW(
        MathLib::PiecewiseLinearInterpolation(
            std::vector<double>{2, 1}, std::vector<double>{20, 10}, true),
        std::runtime_error);
}
