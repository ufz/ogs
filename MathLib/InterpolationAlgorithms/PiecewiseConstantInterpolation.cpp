// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PiecewiseConstantInterpolation.h"

#include <algorithm>
#include <iterator>
#include <range/v3/range/conversion.hpp>
#include <utility>

#include "SupportingPoints.h"

namespace MathLib
{
PiecewiseConstantInterpolation::PiecewiseConstantInterpolation(
    std::vector<double> supporting_points,
    std::vector<double>
        values_at_supp_pnts)
    : supp_pnts_(std::move(supporting_points)),
      values_at_supp_pnts_(std::move(values_at_supp_pnts))
{
    sortAndCheckSupportingPoints(supp_pnts_, values_at_supp_pnts_, 1,
                                 "piecewise constant interpolation", false);
}

PiecewiseConstantInterpolation::PiecewiseConstantInterpolation(
    std::vector<int> const& supporting_points,
    std::vector<double>
        values_at_supp_pnts)
    : PiecewiseConstantInterpolation(
          supporting_points | ranges::to<std::vector<double>>(),
          std::move(values_at_supp_pnts))
{
}

double PiecewiseConstantInterpolation::value(
    double const pnt_to_interpolate) const
{
    if (pnt_to_interpolate <= supp_pnts_.front())
    {
        return values_at_supp_pnts_.front();
    }

    if (supp_pnts_.back() <= pnt_to_interpolate)
    {
        return values_at_supp_pnts_.back();
    }

    auto const& it(std::upper_bound(supp_pnts_.begin(), supp_pnts_.end(),
                                    pnt_to_interpolate));
    // Here access the iterator it without checking is okay since the
    // corner cases are checked above.
    auto const interval_idx = std::distance(supp_pnts_.begin(), it) - 1;

    return values_at_supp_pnts_[interval_idx];
}
}  // end namespace MathLib
