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

    // compute linear interpolation polynom: y = m * (x - support[i]) + value[i]
    const double m((_values_at_supp_pnts[interval_idx + 1] -
                    _values_at_supp_pnts[interval_idx]) /
                   (_supp_pnts[interval_idx + 1] - _supp_pnts[interval_idx]));

    return m * (pnt_to_interpolate - _supp_pnts[interval_idx]) +
           _values_at_supp_pnts[interval_idx];
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
    std::size_t const interval_idx = std::distance(_supp_pnts.begin(), it) - 1;

    // interval_idx = interval_max - 1 - interval_idx;
    if (interval_idx > 1 && interval_idx < _supp_pnts.size() - 2)
    {
        double const tangent_right =
            (_values_at_supp_pnts[interval_idx] -
             _values_at_supp_pnts[interval_idx + 2]) /
            (_supp_pnts[interval_idx] - _supp_pnts[interval_idx + 2]);
        double const tangent_left =
            (_values_at_supp_pnts[interval_idx - 1] -
             _values_at_supp_pnts[interval_idx + 1]) /
            (_supp_pnts[interval_idx - 1] - _supp_pnts[interval_idx + 1]);
        double const w =
            (pnt_to_interpolate - _supp_pnts[interval_idx + 1]) /
            (_supp_pnts[interval_idx] - _supp_pnts[interval_idx + 1]);
        return (1. - w) * tangent_right + w * tangent_left;
    }
    else
    {
        return (_values_at_supp_pnts[interval_idx] -
                _values_at_supp_pnts[interval_idx + 1]) /
               (_supp_pnts[interval_idx] - _supp_pnts[interval_idx + 1]);
    }
}
}  // end MathLib
