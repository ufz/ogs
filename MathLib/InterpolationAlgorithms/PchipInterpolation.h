// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <boost/math/interpolators/pchip.hpp>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace MathLib
{
/**
 * One dimensional shape-preserving piecewise cubic Hermite interpolation
 * (PCHIP), a thin wrapper around `boost::math::interpolators::pchip`.
 *
 * The interpolant passes through all supporting points and is continuously
 * differentiable. The interior nodal slopes are the weighted harmonic means of
 * the neighbouring secant slopes of Fritsch and Butland (1984); they are zero
 * at local extrema and plateaus of the data. Therefore the interpolant is
 * monotone on every interval where the data is monotone, and it does not over-
 * or undershoot the data.
 *
 * At the end points the slope of the adjacent secant is used (the Boost
 * default) unless an end slope is given. A given end slope must have the sign
 * of the adjacent secant and at most three times its magnitude, which keeps
 * the end interval monotone; otherwise it is a fatal error.
 *
 * At least four supporting points are required (Boost requirement).
 *
 * As in PiecewiseLinearInterpolation, the value is constant and the
 * derivative is zero outside of the range of the supporting points.
 */
class PchipInterpolation final
{
public:
    /// The supporting points need not be sorted but must be distinct. Fewer
    /// than four supporting points are a fatal error.
    /// \param curve_name Name of the curve, used in error messages.
    /// \param supporting_points Abscissae of the supporting points.
    /// \param values_at_supp_pnts Values at the supporting points.
    /// \param left_end_slope Slope of the interpolant at the smallest
    /// supporting point; the secant slope of the first interval if not given.
    /// It must lie between 0 and three times that secant slope.
    /// \param right_end_slope Slope of the interpolant at the largest
    /// supporting point; the secant slope of the last interval if not given.
    /// It must lie between 0 and three times that secant slope.
    PchipInterpolation(std::string_view curve_name,
                       std::vector<double>
                           supporting_points,
                       std::vector<double>
                           values_at_supp_pnts,
                       std::optional<double>
                           left_end_slope,
                       std::optional<double>
                           right_end_slope);

    /// Evaluates the cubic Hermite polynomial. Points outside of
    /// \f$[x_{\min}, x_{\max}]\f$ are set to \f$x_{\min}\f$ or \f$x_{\max}\f$.
    double getValue(double pnt_to_interpolate) const;

    /// Derivative of the cubic Hermite polynomial; 0 outside of
    /// \f$[x_{\min}, x_{\max}]\f$.
    double getDerivative(double pnt_to_interpolate) const;

    double getSupportMin() const;
    double getSupportMax() const;

private:
    /// \param sorted_points Supporting points sorted in ascending order and
    /// the values at them, as returned by sortAndCheckSupportingPoints().
    PchipInterpolation(std::string_view curve_name,
                       std::pair<std::vector<double>, std::vector<double>>
                           sorted_points,
                       std::optional<double>
                           left_end_slope,
                       std::optional<double>
                           right_end_slope);

    // Stored here since boost's pchip provides no domain() before 1.92.
    double x_min_;
    double x_max_;
    boost::math::interpolators::pchip<std::vector<double>> pchip_;
};
}  // namespace MathLib
