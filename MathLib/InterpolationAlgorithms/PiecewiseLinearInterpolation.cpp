/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-07
 * \brief  Implementation of the PiecewiseLinearInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "PiecewiseLinearInterpolation.h"

#include <cmath>
#include <limits>
#include <utility>

#include "BaseLib/Error.h"
#include "BaseLib/quicksort.h"

namespace MathLib
{
PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(
    std::vector<double> supporting_points,
    std::vector<double>
        values_at_supp_pnts,
    bool supp_pnts_sorted)
    : supp_pnts_(std::move(supporting_points)),
      values_at_supp_pnts_(std::move(values_at_supp_pnts))
{
    if (!supp_pnts_sorted)
    {
        BaseLib::quicksort<double, double>(
            supp_pnts_, static_cast<std::size_t>(0), supp_pnts_.size(),
            values_at_supp_pnts_);
    }

    const auto it = std::adjacent_find(supp_pnts_.begin(), supp_pnts_.end());
    if (it != supp_pnts_.end())
    {
        const std::size_t i = std::distance(supp_pnts_.begin(), it);
        OGS_FATAL(
            "Variable {:d} and variable {:d} are the same. Piecewise linear "
            "interpolation is not possible\n",
            i, i + 1);
    }
}

double PiecewiseLinearInterpolation::getValue(double pnt_to_interpolate) const
{
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
    const double m = (f_r - f) / (x_r - x);

    return m * (pnt_to_interpolate - x) + f;
}

double PiecewiseLinearInterpolation::getDerivative(
    double const pnt_to_interpolate) const
{
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
