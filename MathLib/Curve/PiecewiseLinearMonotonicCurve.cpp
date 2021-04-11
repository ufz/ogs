/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
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
    std::vector<double> x, std::vector<double> y)
    : PiecewiseLinearInterpolation(std::move(x), std::move(y), false)
{
    if (!isStrongMonotonic())
    {
        OGS_FATAL("The given curve is not strong monotonic.");
    }
}

bool PiecewiseLinearMonotonicCurve::isStrongMonotonic() const
{
    const double gradient0 = getDerivative(supp_pnts_[0]);

    if (std::abs(gradient0) < std::numeric_limits<double>::min())
    {
        return false;
    }

    return std::none_of(supp_pnts_.begin(), supp_pnts_.end(),
                        [&](const double p) {
                            return this->getDerivative(p) * gradient0 <= 0.;
                        });
}

double PiecewiseLinearMonotonicCurve::getInverseVariable(const double y) const
{
    std::size_t interval_idx = 0;
    if (values_at_supp_pnts_.front() < values_at_supp_pnts_.back())
    {
        if (y <= values_at_supp_pnts_.front())
        {
            return supp_pnts_[0];
        }
        if (y >= values_at_supp_pnts_.back())
        {
            return supp_pnts_[supp_pnts_.size() - 1];
        }

        // search interval that has the point inside
        auto const& it(std::lower_bound(values_at_supp_pnts_.begin(),
                                        values_at_supp_pnts_.end(), y));
        interval_idx = std::distance(values_at_supp_pnts_.begin(), it) - 1;
    }
    else
    {
        if (y >= values_at_supp_pnts_.front())
        {
            return supp_pnts_[0];
        }
        if (y <= values_at_supp_pnts_.back())
        {
            return supp_pnts_[supp_pnts_.size() - 1];
        }

        // search interval in the reverse direction for the point inside
        auto const& it(std::lower_bound(values_at_supp_pnts_.rbegin(),
                                        values_at_supp_pnts_.rend(), y));
        interval_idx = values_at_supp_pnts_.size() -
                       (std::distance(values_at_supp_pnts_.rbegin(), it)) - 1;
    }

    const double xi_1 = supp_pnts_[interval_idx + 1];
    const double xi = supp_pnts_[interval_idx];
    const double yi_1 = values_at_supp_pnts_[interval_idx + 1];
    const double yi = values_at_supp_pnts_[interval_idx];

    // compute gradient: m = (x_{i+1} - x_i) / (y_{i+1} - y_i)
    const double m = (xi_1 - xi) / (yi_1 - yi);

    // compute the variable by linear interpolation:  x = m * (y - y_i) + x_i,
    // and then return the result.
    return m * (y - yi) + xi;
}

}  // namespace MathLib
