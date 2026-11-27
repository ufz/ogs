// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PiecewiseLinearInterpolation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <range/v3/range/conversion.hpp>
#include <utility>

#include "SupportingPoints.h"

namespace MathLib
{
PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(
    std::vector<double> supporting_points,
    std::vector<double>
        values_at_supp_pnts,
    bool supp_pnts_sorted,
    std::string_view const description)
    : supp_pnts_(std::move(supporting_points)),
      values_at_supp_pnts_(std::move(values_at_supp_pnts))
{
    sortAndCheckSupportingPoints(supp_pnts_, values_at_supp_pnts_, 1,
                                 description, supp_pnts_sorted);
}

PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(
    std::vector<int> const& supporting_points,
    std::vector<double> const& values_at_supp_pnts,
    bool supp_pnts_sorted)
    : supp_pnts_(supporting_points | ranges::to<std::vector<double>>()),
      values_at_supp_pnts_(values_at_supp_pnts)
{
    sortAndCheckSupportingPoints(supp_pnts_, values_at_supp_pnts_, 1,
                                 "piecewise linear interpolation",
                                 supp_pnts_sorted);
}

double PiecewiseLinearInterpolation::getValue(double pnt_to_interpolate) const
{
    // All comparisons are false for NaN, so the search below would give the
    // interval index -1.
    if (std::isnan(pnt_to_interpolate))
    {
        return pnt_to_interpolate;
    }

    // search interval that has the point inside
    if (pnt_to_interpolate <= supp_pnts_.front())
    {
        return values_at_supp_pnts_[0];
    }

    if (supp_pnts_.back() <= pnt_to_interpolate)
    {
        return values_at_supp_pnts_[supp_pnts_.size() - 1];
    }

    auto const& it(std::lower_bound(supp_pnts_.begin(), supp_pnts_.end(),
                                    pnt_to_interpolate));
    std::size_t const interval_idx = std::distance(supp_pnts_.begin(), it) - 1;

    // support points.
    double const x = supp_pnts_[interval_idx];
    double const x_r = supp_pnts_[interval_idx + 1];

    // values.
    double const f = values_at_supp_pnts_[interval_idx];
    double const f_r = values_at_supp_pnts_[interval_idx + 1];

    // compute linear interpolation polynom: y = m * (x - support[i]) + value[i]
    const double t = (pnt_to_interpolate - x) / (x_r - x);
    return std::lerp(f, f_r, t);
}

double PiecewiseLinearInterpolation::getDerivative(
    double const pnt_to_interpolate) const
{
    // All comparisons are false for NaN, so the search below would give the
    // interval index 0 and read before the first value.
    if (std::isnan(pnt_to_interpolate))
    {
        return pnt_to_interpolate;
    }

    // A single support point describes a constant function.
    if (supp_pnts_.size() < 2)
    {
        return 0;
    }

    if (pnt_to_interpolate < supp_pnts_.front() ||
        supp_pnts_.back() < pnt_to_interpolate)
    {
        return 0;
    }

    auto const& it(std::lower_bound(supp_pnts_.begin(), supp_pnts_.end(),
                                    pnt_to_interpolate));
    std::size_t interval_idx = std::distance(supp_pnts_.begin(), it);

    if (pnt_to_interpolate == supp_pnts_.front())
    {
        interval_idx = 1;
    }

    if (interval_idx > 1 && interval_idx < supp_pnts_.size() - 2)
    {
        // left and right support points.
        double const x_ll = supp_pnts_[interval_idx - 2];
        double const x_l = supp_pnts_[interval_idx - 1];
        double const x = supp_pnts_[interval_idx];
        double const x_r = supp_pnts_[interval_idx + 1];

        // left and right values.
        double const f_ll = values_at_supp_pnts_[interval_idx - 2];
        double const f_l = values_at_supp_pnts_[interval_idx - 1];
        double const f = values_at_supp_pnts_[interval_idx];
        double const f_r = values_at_supp_pnts_[interval_idx + 1];

        double const tangent_right = (f_l - f_r) / (x_l - x_r);
        double const tangent_left = (f_ll - f) / (x_ll - x);
        double const w = (pnt_to_interpolate - x) / (x_l - x);
        return (1. - w) * tangent_right + w * tangent_left;
    }

    return (values_at_supp_pnts_[interval_idx] -
            values_at_supp_pnts_[interval_idx - 1]) /
           (supp_pnts_[interval_idx] - supp_pnts_[interval_idx - 1]);
}

double PiecewiseLinearInterpolation::getSupportMax() const
{
    assert(!supp_pnts_.empty());
    return supp_pnts_.back();
}
double PiecewiseLinearInterpolation::getSupportMin() const
{
    assert(!supp_pnts_.empty());
    return supp_pnts_.front();
}
}  // namespace MathLib
