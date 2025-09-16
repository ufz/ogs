/**
 * \file
 * \brief  Definition of the PiecewiseConstantInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/Logging.h"

namespace MathLib
{
/**
 * This class implements a one dimensional piecewise constant interpolation
 * algorithm.
 */
template <typename T = double>
class PiecewiseConstantInterpolation
{
public:
    /**
     * The constructor stores the vector of supporting points
     * \f$(x_0, x_1, \dots, x_n)\f$ and the entries of the vector of values at
     * the supporting points \f$(y_0, y_1, \dots, y_n)\f$ where \f$n\f$
     * is the number of entries of the vector. The number of supporting
     * points must be equal to the number of values at the supporting points.
     * It is assumed that \f$x_j\f$ corresponds to \f$y_j\f$ for all \f$j \in
     * [0, n]\f$.
     *
     * Furthermore, it is assumed that the supporting points are sorted, i.e.
     * \f$x_0 < x_1 < \dots < x_n\f$.
     * @param supporting_points vector of supporting points
     * @param values_at_supp_pnts vector of values at the supporting points
     * one can set the switch to true
     */
    PiecewiseConstantInterpolation(
        std::vector<T> const& supporting_points,
        std::vector<double> const& values_at_supp_pnts)
        : supp_pnts_(supporting_points),
          values_at_supp_pnts_(values_at_supp_pnts)
    {
        if (supp_pnts_.size() != values_at_supp_pnts_.size())
        {
            OGS_FATAL(
                "Inconsistent data given to PiecewiseConstantInterpolation, "
                "number of given supporting points is {}, number of given "
                "values is {}.",
                supp_pnts_.size(), values_at_supp_pnts_.size());
        }
        if (supp_pnts_.empty())
        {
            ERR("PiecewiseConstantInterpolation: passed empty vector.");
        }
    }

    /**
     * \brief Calculates the interpolation value.
     * @param pnt_to_interpolate The point the interpolation value is calculated
     * for.
     * If the `pnt_to_interpolate` is outside the interval \f$[x_{\min},
     * x_{\max}]\f$, where \f$x_{\min} = \min_{1 \le j \le n} x_j\f$ and
     * \f$x_{\max} = \max_{1 \le j \le n} x_j\f$. If `point_to_interpolate` is
     * smaller than \f$x_{\min}\f$ then \f$y_{\min}\f$ is returned. Analogously,
     * if the `point_to_interpolate` is greater than \f$x_{\max}\f$ then
     * \f$y_{\max}\f$ is returned.
     * @return The interpolated value.
     */
    double value(double const pnt_to_interpolate) const
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

private:
    std::vector<T> const& supp_pnts_;
    std::vector<double> const& values_at_supp_pnts_;
};
}  // end namespace MathLib
