/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-07
 * \brief  Implementation of the PiecewiseLinearInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <cmath>
#include <limits>
#include <utility>

#include "BaseLib/quicksort.h"
#include "PiecewiseLinearInterpolation.h"

namespace MathLib
{
PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(
    std::vector<double>&& supporting_points,
    std::vector<double>&& values_at_supp_pnts,
    bool supp_pnts_sorted)
    : _supp_pnts(std::move(supporting_points)),
      _values_at_supp_pnts(std::move(values_at_supp_pnts))
{
    if (!supp_pnts_sorted)
    {
        BaseLib::quicksort<double, double>(
            _supp_pnts, static_cast<std::size_t>(0), _supp_pnts.size(),
            _values_at_supp_pnts);
    }
}

double PiecewiseLinearInterpolation::getValue(double pnt_to_interpolate) const
{
    // search interval that has the point inside
    if (pnt_to_interpolate <= _supp_pnts.front())
    {
        return _values_at_supp_pnts[0];
    }

    if (_supp_pnts.back() <= pnt_to_interpolate)
    {
        return _values_at_supp_pnts[_supp_pnts.size() - 1];
    }

    auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(),
                                    pnt_to_interpolate));
    std::size_t const interval_idx = std::distance(_supp_pnts.begin(), it) - 1;

    // support points.
    double const x = _supp_pnts[interval_idx];
    double const x_r = _supp_pnts[interval_idx + 1];

    // values.
    double const f = _values_at_supp_pnts[interval_idx];
    double const f_r = _values_at_supp_pnts[interval_idx + 1];

    // compute linear interpolation polynom: y = m * (x - support[i]) + value[i]
    const double m = (f_r - f) / (x_r - x);

    return m * (pnt_to_interpolate - x) + f;
}

double PiecewiseLinearInterpolation::getDerivative(
    double const pnt_to_interpolate) const
{
    if (pnt_to_interpolate < _supp_pnts.front() ||
        _supp_pnts.back() < pnt_to_interpolate)
    {
        return 0;
    }

    auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(),
                                    pnt_to_interpolate));
    std::size_t interval_idx = std::distance(_supp_pnts.begin(), it);

    if (pnt_to_interpolate == _supp_pnts.front())
    {
        interval_idx = 1;
    }

    if (interval_idx > 2 && interval_idx < _supp_pnts.size() - 1)
    {
        // left and right support points.
        double const x_ll = _supp_pnts[interval_idx - 2];
        double const x_l = _supp_pnts[interval_idx - 1];
        double const x = _supp_pnts[interval_idx];
        double const x_r = _supp_pnts[interval_idx + 1];

        // left and right values.
        double const f_ll = _values_at_supp_pnts[interval_idx - 2];
        double const f_l = _values_at_supp_pnts[interval_idx - 1];
        double const f = _values_at_supp_pnts[interval_idx];
        double const f_r = _values_at_supp_pnts[interval_idx + 1];

        double const tangent_right = (f_l - f_r) / (x_l - x_r);
        double const tangent_left = (f_ll - f) / (x_ll - x);
        double const w = (pnt_to_interpolate - x) / (x_l - x);
        return (1. - w) * tangent_right + w * tangent_left;
    }
    else
    {
        return (_values_at_supp_pnts[interval_idx] -
                _values_at_supp_pnts[interval_idx - 1]) /
               (_supp_pnts[interval_idx] - _supp_pnts[interval_idx - 1]);
    }
}

double PiecewiseLinearInterpolation::getSupportMax() const
{
    assert(!_supp_pnts.empty());
    return _supp_pnts.back();
}
double PiecewiseLinearInterpolation::getSupportMin() const
{
    assert(!_supp_pnts.empty());
    return _supp_pnts.front();
}
}  // end MathLib
