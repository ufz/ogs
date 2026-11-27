// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PchipInterpolation.h"

namespace
{
// Data and reference values computed with
// scipy.interpolate.CubicHermiteSpline(x, y, slopes), where the slopes are the
// ones of boost::math::interpolators::pchip: Fritsch-Butland harmonic means in
// the interior (as scipy.interpolate.PchipInterpolator) and the secant slope
// of the adjacent interval at the two end points.
std::vector<double> const supporting_points = {0, 1, 2.5, 3, 5, 6.5};
std::vector<double> const values = {0, 1, 1.2, 4, 4.1, 2.0};
}  // namespace

class MathLibPchipInterpolation : public ::testing::Test
{
protected:
    // End slopes default to the secant slopes of the end intervals.
    MathLib::PchipInterpolation const default_pchip{
        "test", supporting_points, values, std::nullopt, std::nullopt};
};

TEST_F(MathLibPchipInterpolation, PassesThroughSupportingPoints)
{
    for (std::size_t i = 0; i < supporting_points.size(); ++i)
    {
        EXPECT_NEAR(values[i], default_pchip.getValue(supporting_points[i]),
                    1e-14);
    }
}

TEST_F(MathLibPchipInterpolation, ReferenceValuesAndDerivatives)
{
    std::vector<double> const q = {0.5, 1.7, 2.75, 4, 5.5, 6.2};
    std::vector<double> const value = {0.59400826446281, 1.08542851782734,
                                       2.61164558760836, 4.08083700440529,
                                       3.71111111111111, 2.4872};
    std::vector<double> const derivative = {1.18801652892562,
                                            0.0595088717082851,
                                            8.291743640756,
                                            0.0441629955947135,
                                            -1.4,
                                            -1.792};
    for (std::size_t i = 0; i < q.size(); ++i)
    {
        EXPECT_NEAR(value[i], default_pchip.getValue(q[i]), 1e-12) << q[i];
        EXPECT_NEAR(derivative[i], default_pchip.getDerivative(q[i]), 1e-12)
            << q[i];
    }
}

TEST_F(MathLibPchipInterpolation, GivenEndSlopes)
{
    // Slope 0 at the first and -0.5 at the last point, otherwise the slopes
    // of boost::math::interpolators::pchip (reference values from a cubic
    // Hermite spline with these slopes).
    MathLib::PchipInterpolation const pchip("test", supporting_points, values,
                                            0.0, -0.5);
    EXPECT_NEAR(0.0, pchip.getDerivative(supporting_points.front()), 1e-14);
    EXPECT_NEAR(-0.5, pchip.getDerivative(supporting_points.back()), 1e-14);

    std::vector<double> const q = {0.5, 1.7, 4, 5.5, 6.2};
    std::vector<double> const value = {0.46900826446281, 1.08542851782734,
                                       4.08083700440529, 3.61111111111111,
                                       2.3144};
    std::vector<double> const derivative = {
        1.43801652892562, 0.0595088717082851, 0.0441629955947135, -1.7, -1.504};
    for (std::size_t i = 0; i < q.size(); ++i)
    {
        EXPECT_NEAR(value[i], pchip.getValue(q[i]), 1e-12) << q[i];
        EXPECT_NEAR(derivative[i], pchip.getDerivative(q[i]), 1e-12) << q[i];
    }

    // Without end slopes the default secant slope is used.
    EXPECT_NEAR(1.0, default_pchip.getDerivative(supporting_points.front()),
                1e-14);
}

TEST_F(MathLibPchipInterpolation,
       ZeroSlopeAtFirstPointGivesContinuousDerivative)
{
    // Constant extrapolation below the first point has slope 0; with a zero
    // end slope the derivative is continuous there.
    MathLib::PchipInterpolation const pchip("test", supporting_points, values,
                                            0.0, std::nullopt);
    double const eps = 1e-9;
    EXPECT_NEAR(pchip.getDerivative(supporting_points.front() - eps),
                pchip.getDerivative(supporting_points.front() + eps), 1e-8);
}

TEST_F(MathLibPchipInterpolation, DerivativeIsDerivativeOfValue)
{
    double const eps = 1e-6;
    for (double q = 0.01; q < 6.5; q += 0.0731)
    {
        double const fd = (default_pchip.getValue(q + eps) -
                           default_pchip.getValue(q - eps)) /
                          (2 * eps);
        EXPECT_NEAR(fd, default_pchip.getDerivative(q), 1e-6) << q;
    }
}

TEST_F(MathLibPchipInterpolation, ContinuousDerivativeAtSupportingPoints)
{
    double const eps = 1e-10;
    for (std::size_t i = 1; i + 1 < supporting_points.size(); ++i)
    {
        EXPECT_NEAR(default_pchip.getDerivative(supporting_points[i] - eps),
                    default_pchip.getDerivative(supporting_points[i] + eps),
                    1e-8)
            << supporting_points[i];
    }
}

TEST_F(MathLibPchipInterpolation, PreservesMonotonicity)
{
    // Steep increase followed by a plateau, which makes a cubic spline
    // overshoot.
    MathLib::PchipInterpolation const pchip("test", {0, 1, 2, 3, 4, 5},
                                            {0, 0.01, 1, 1, 1.01, 3},
                                            std::nullopt, std::nullopt);
    double previous = pchip.getValue(0);
    for (double q = 0; q <= 5; q += 0.001)
    {
        double const value = pchip.getValue(q);
        EXPECT_GE(value, previous - 1e-14) << q;
        EXPECT_GE(pchip.getDerivative(q), -1e-12) << q;
        EXPECT_GE(value, 0.);
        EXPECT_LE(value, 3.);
        previous = value;
    }
    // Plateau: flat between the equal values.
    EXPECT_NEAR(1.0, pchip.getValue(2.5), 1e-14);
    EXPECT_NEAR(0.0, pchip.getDerivative(2.5), 1e-14);
}

TEST_F(MathLibPchipInterpolation,
       PreservesMonotonicityOfDecreasingDataWithZeroLeftEndSlope)
{
    // Retention curve of the RichardsFlow_2d_small benchmark.
    std::vector<double> const capillary_pressure = {
        1421.99634, 1707.55369, 2050.45507, 2462.21599, 2892.29884, 3232.80969,
        3590.71602, 3988.24638, 4429.78757, 4920.21206, 5464.93174, 9238.24221,
        10261.0136, 12658.7874, 15616.8683, 26399.6732, 44627.5611};
    std::vector<double> const saturation = {
        1,          0.98684211, 0.97368421, 0.96052632, 0.94894737, 0.92105263,
        0.89473684, 0.86842105, 0.84210526, 0.81578947, 0.78947368, 0.65789474,
        0.63157895, 0.57894737, 0.52631579, 0.39473684, 0.26315789};
    MathLib::PchipInterpolation const pchip("test", capillary_pressure,
                                            saturation, 0.0, std::nullopt);
    EXPECT_NEAR(0.0, pchip.getDerivative(capillary_pressure.front()), 1e-14);

    double previous = pchip.getValue(capillary_pressure.front());
    for (double q = capillary_pressure.front(); q <= capillary_pressure.back();
         q += 10)
    {
        double const value = pchip.getValue(q);
        EXPECT_LE(value, previous + 1e-14) << q;
        EXPECT_LE(pchip.getDerivative(q), 1e-15) << q;
        EXPECT_LE(value, saturation.front() + 1e-14) << q;
        EXPECT_GE(value, saturation.back() - 1e-14) << q;
        previous = value;
    }
}

TEST_F(MathLibPchipInterpolation, ReproducesLinearData)
{
    MathLib::PchipInterpolation const pchip("test", {-1, 0, 0.5, 3, 4},
                                            {-3, -1, 0, 5, 7}, std::nullopt,
                                            std::nullopt);
    for (double q = -1; q <= 4; q += 0.0937)
    {
        EXPECT_NEAR(2 * q - 1, pchip.getValue(q), 1e-13) << q;
        EXPECT_NEAR(2., pchip.getDerivative(q), 1e-13) << q;
    }
}

TEST_F(MathLibPchipInterpolation, ConstantOutsideOfSupport)
{
    EXPECT_EQ(values.front(), default_pchip.getValue(-10));
    EXPECT_EQ(values.back(), default_pchip.getValue(10));
    EXPECT_EQ(0., default_pchip.getDerivative(-10));
    EXPECT_EQ(0., default_pchip.getDerivative(10));
    EXPECT_EQ(supporting_points.front(), default_pchip.getSupportMin());
    EXPECT_EQ(supporting_points.back(), default_pchip.getSupportMax());
}

TEST_F(MathLibPchipInterpolation, NanArgumentGivesNan)
{
    // Boost reads past the end of its data for a NaN argument (reported as
    // heap-buffer-overflow by AddressSanitizer); the result is NaN anyway.
    EXPECT_TRUE(std::isnan(
        default_pchip.getValue(std::numeric_limits<double>::quiet_NaN())));
    EXPECT_TRUE(std::isnan(
        default_pchip.getDerivative(std::numeric_limits<double>::quiet_NaN())));
}

TEST_F(MathLibPchipInterpolation, UnsortedSupportingPoints)
{
    std::vector<double> const xr(supporting_points.rbegin(),
                                 supporting_points.rend());
    std::vector<double> const yr(values.rbegin(), values.rend());
    MathLib::PchipInterpolation const reversed("test", xr, yr, std::nullopt,
                                               std::nullopt);
    for (double q = 0; q <= 6.5; q += 0.13)
    {
        EXPECT_EQ(default_pchip.getValue(q), reversed.getValue(q));
        EXPECT_EQ(default_pchip.getDerivative(q), reversed.getDerivative(q));
    }
}

TEST_F(MathLibPchipInterpolation, FewerThanFourPointsAreRejected)
{
    EXPECT_THROW(MathLib::PchipInterpolation("test", {1, 3}, {2, 6},
                                             std::nullopt, std::nullopt),
                 std::runtime_error);
    EXPECT_THROW(MathLib::PchipInterpolation("test", {1, 2, 3}, {2, 6, 3},
                                             std::nullopt, std::nullopt),
                 std::runtime_error);
    EXPECT_NO_THROW(MathLib::PchipInterpolation(
        "test", {1, 2, 3, 4}, {2, 6, 3, 4}, std::nullopt, std::nullopt));
}

TEST_F(MathLibPchipInterpolation, DuplicateSupportingPointsAreRejected)
{
    EXPECT_THROW(MathLib::PchipInterpolation("test", {1, 2, 2, 4}, {2, 6, 3, 4},
                                             std::nullopt, std::nullopt),
                 std::runtime_error);
}

TEST_F(MathLibPchipInterpolation, EndSlopesBreakingMonotonicityAreRejected)
{
    // Secant slopes of the end intervals: 1 on the left, -1.4 on the right.
    EXPECT_THROW(MathLib::PchipInterpolation("test", supporting_points, values,
                                             -0.1, std::nullopt),
                 std::runtime_error);
    EXPECT_THROW(MathLib::PchipInterpolation("test", supporting_points, values,
                                             3.1, std::nullopt),
                 std::runtime_error);
    EXPECT_THROW(MathLib::PchipInterpolation("test", supporting_points, values,
                                             std::nullopt, 0.1),
                 std::runtime_error);
    EXPECT_THROW(MathLib::PchipInterpolation("test", supporting_points, values,
                                             std::nullopt, -4.3),
                 std::runtime_error);
    EXPECT_NO_THROW(MathLib::PchipInterpolation("test", supporting_points,
                                                values, 0.0, 0.0));
    EXPECT_NO_THROW(MathLib::PchipInterpolation("test", supporting_points,
                                                values, 3.0, -4.1));
    // Exactly three times the secant slopes despite round-off in the computed
    // right secant (-1.3999999999999997).
    EXPECT_NO_THROW(MathLib::PchipInterpolation("test", supporting_points,
                                                values, 3.0, -4.2));

    // Flat first interval: only a zero slope keeps it flat.
    std::vector<double> const flat = {1, 1, 2, 3};
    std::vector<double> const xs = {0, 1, 2, 3};
    EXPECT_THROW(
        MathLib::PchipInterpolation("test", xs, flat, 0.1, std::nullopt),
        std::runtime_error);
    EXPECT_NO_THROW(
        MathLib::PchipInterpolation("test", xs, flat, 0.0, std::nullopt));
}

TEST_F(MathLibPchipInterpolation, LimitingEndSlopesKeepEndIntervalsWithinData)
{
    MathLib::PchipInterpolation const pchip("test", supporting_points, values,
                                            3.0, -4.1);
    for (double q = supporting_points[0]; q <= supporting_points[1]; q += 0.001)
    {
        EXPECT_GE(pchip.getValue(q), values[0] - 1e-14) << q;
        EXPECT_LE(pchip.getValue(q), values[1] + 1e-14) << q;
    }
    for (double q = supporting_points[4]; q <= supporting_points[5]; q += 0.001)
    {
        EXPECT_GE(pchip.getValue(q), values[5] - 1e-14) << q;
        EXPECT_LE(pchip.getValue(q), values[4] + 1e-14) << q;
    }
}
