// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PchipInterpolation.h"

#include <cmath>
#include <limits>
#include <string_view>
#include <utility>

#include "BaseLib/Error.h"
#include "SupportingPoints.h"

namespace MathLib
{
namespace
{
void checkEndSlope(std::string_view const curve_name,
                   std::optional<double> const slope, double const secant,
                   std::string_view const end)
{
    if (!slope)
    {
        return;
    }
    // The interior slopes are at most three times the adjacent secant slopes.
    // An end slope of the sign of the secant and at most three times its
    // magnitude therefore keeps the cubic on the end interval monotone, which
    // is the sufficient condition of Fritsch and Carlson (1980). The upper
    // bound allows for a few ulp of round-off in the computed secant, so that
    // exactly three times the secant of the input data is accepted.
    double const max_ratio =
        3 * (1 + 4 * std::numeric_limits<double>::epsilon());
    bool const is_monotone =
        secant == 0 ? *slope == 0
                    : (0 <= *slope / secant && *slope / secant <= max_ratio);
    if (!is_monotone)
    {
        OGS_FATAL(
            "The {:s} end slope {:g} of the PCHIP curve '{:s}' must have the "
            "sign of the slope {:g} of the adjacent interval and at most three "
            "times its magnitude (0 if the interval is flat). Otherwise the "
            "interpolant over- or undershoots the supporting values there.",
            end, *slope, curve_name, secant);
    }
}

std::pair<std::vector<double>, std::vector<double>> sortAndCheck(
    std::string_view const curve_name, std::vector<double> x,
    std::vector<double> y)
{
    // Four supporting points are required by boost's pchip.
    sortAndCheckSupportingPoints(
        x, y, 4, fmt::format("PCHIP curve '{:s}'", curve_name), false);
    return {std::move(x), std::move(y)};
}

/// \pre x is sorted and has at least four entries, see sortAndCheck().
boost::math::interpolators::pchip<std::vector<double>> createPchip(
    std::string_view const curve_name, std::vector<double> x,
    std::vector<double> y, std::optional<double> const left_end_slope,
    std::optional<double> const right_end_slope)
{
    auto const n = x.size();
    checkEndSlope(curve_name, left_end_slope, (y[1] - y[0]) / (x[1] - x[0]),
                  "left");
    checkEndSlope(curve_name, right_end_slope,
                  (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]), "right");

    // Boost takes NaN for "use the secant slope".
    return boost::math::interpolators::pchip<std::vector<double>>(
        std::move(x), std::move(y),
        left_end_slope.value_or(std::numeric_limits<double>::quiet_NaN()),
        right_end_slope.value_or(std::numeric_limits<double>::quiet_NaN()));
}
}  // namespace

PchipInterpolation::PchipInterpolation(
    std::string_view const curve_name,
    std::vector<double>
        supporting_points,
    std::vector<double>
        values_at_supp_pnts,
    std::optional<double> const left_end_slope,
    std::optional<double> const right_end_slope)
    : PchipInterpolation(curve_name,
                         sortAndCheck(curve_name, std::move(supporting_points),
                                      std::move(values_at_supp_pnts)),
                         left_end_slope, right_end_slope)
{
}

PchipInterpolation::PchipInterpolation(
    std::string_view const curve_name,
    std::pair<std::vector<double>, std::vector<double>>
        sorted_points,
    std::optional<double> const left_end_slope,
    std::optional<double> const right_end_slope)
    // x_min_ and x_max_ are initialized before the points are moved into
    // pchip_, see the order of the members.
    : x_min_(sorted_points.first.front()),
      x_max_(sorted_points.first.back()),
      pchip_(createPchip(curve_name, std::move(sorted_points.first),
                         std::move(sorted_points.second), left_end_slope,
                         right_end_slope))
{
}

double PchipInterpolation::getValue(double const pnt_to_interpolate) const
{
    // Boost reads past its data for NaN, since all comparisons are false.
    if (std::isnan(pnt_to_interpolate))
    {
        return pnt_to_interpolate;
    }
    if (pnt_to_interpolate <= x_min_)
    {
        return pchip_(x_min_);
    }
    if (x_max_ <= pnt_to_interpolate)
    {
        return pchip_(x_max_);
    }
    return pchip_(pnt_to_interpolate);
}

double PchipInterpolation::getDerivative(double const pnt_to_interpolate) const
{
    if (std::isnan(pnt_to_interpolate))
    {
        return pnt_to_interpolate;
    }
    if (pnt_to_interpolate < x_min_ || x_max_ < pnt_to_interpolate)
    {
        return 0.;
    }
    return pchip_.prime(pnt_to_interpolate);
}

double PchipInterpolation::getSupportMin() const
{
    return x_min_;
}

double PchipInterpolation::getSupportMax() const
{
    return x_max_;
}
}  // namespace MathLib
