/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   PiecewiseLinearMonotonicCurve.cpp
 *
 * Created on November 11, 2016, 10:49 AM
 */

#include "PiecewiseLinearMonotonicCurve.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "BaseLib/Error.h"

namespace MathLib
{
PiecewiseLinearMonotonicCurve::PiecewiseLinearMonotonicCurve(
    std::vector<double>&& x, std::vector<double>&& y)
    : PiecewiseLinearInterpolation(std::move(x), std::move(y), false)
{
    if (!isStrongMonotonic())
    {
        OGS_FATAL("The given curve is not strong monotonic.");
    }
}

bool PiecewiseLinearMonotonicCurve::isStrongMonotonic() const
{
    const double gradient0 = getDerivative(_supp_pnts[0]);

    if (std::fabs(gradient0) < std::numeric_limits<double>::min())
        return false;
    else
    {
        return std::none_of(_supp_pnts.begin(), _supp_pnts.end(),
                            [&](const double p) {
                                return this->getDerivative(p) * gradient0 <= 0.;
                            });
    }
}

double PiecewiseLinearMonotonicCurve::getInversVariable(const double y) const
{
    std::size_t interval_idx = 0;
    if (_values_at_supp_pnts.front() < _values_at_supp_pnts.back())
    {
        if (y <= _values_at_supp_pnts.front())
        {
            return _supp_pnts[0];
        }
        else if (y >= _values_at_supp_pnts.back())
        {
            return _supp_pnts[_supp_pnts.size() - 1];
        }
        else
        {
            // search interval that has the point inside
            auto const& it(std::lower_bound(_values_at_supp_pnts.begin(),
                                            _values_at_supp_pnts.end(), y));
            interval_idx = std::distance(_values_at_supp_pnts.begin(), it) - 1;
        }
    }
    else
    {
        if (y >= _values_at_supp_pnts.front())
        {
            return _supp_pnts[0];
        }
        else if (y <= _values_at_supp_pnts.back())
        {
            return _supp_pnts[_supp_pnts.size() - 1];
        }
        else
        {
            // search interval in the reverse direction for the point inside
            auto const& it(std::lower_bound(_values_at_supp_pnts.rbegin(),
                                            _values_at_supp_pnts.rend(), y));
            interval_idx = _values_at_supp_pnts.size() -
                           (std::distance(_values_at_supp_pnts.rbegin(), it)) -
                           1;
        }
    }

    const double xi_1 = _supp_pnts[interval_idx + 1];
    const double xi = _supp_pnts[interval_idx];
    const double yi_1 = _values_at_supp_pnts[interval_idx + 1];
    const double yi = _values_at_supp_pnts[interval_idx];

    // compute gradient: m = (x_{i+1} - x_i) / (y_{i+1} - y_i)
    const double m = (xi_1 - xi) / (yi_1 - yi);

    // compute the variable by linear interpolation:  x = m * (y - y_i) + x_i,
    // and then return the result.
    return m * (y - yi) + xi;
}

}  // namespace
