/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   PiecewiseLinearCurve.cpp
 *
 * Created on November 11, 2016, 10:49 AM
 */

#include "PiecewiseLinearCurve.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace MathLib
{
bool PiecewiseLinearCurve::isStrongMonotonic() const
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

double PiecewiseLinearCurve::getVariable(const double y) const
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

    // compute gradient: m = (x_{i+1} - x_i) / (y_{i+1} - y_i)
    const double m = (_supp_pnts[interval_idx + 1] - _supp_pnts[interval_idx]) /
                     (_values_at_supp_pnts[interval_idx + 1] -
                      _values_at_supp_pnts[interval_idx]);

    // compute the variable by linear interpolation:  x = m * (y - y_i) + x_i,
    // and then return the result.
    return m * (y - _values_at_supp_pnts[interval_idx]) +
           _supp_pnts[interval_idx];
}

}  // namespace
