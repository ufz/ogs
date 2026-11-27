// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace MathLib
{
/**
 * This class implements a one dimensional piecewise constant interpolation
 * algorithm.
 */
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
     * The supporting points need not be sorted but must be pairwise different
     * and at least one must be given; otherwise it is a fatal error, see
     * sortAndCheckSupportingPoints().
     * @param supporting_points vector of supporting points
     * @param values_at_supp_pnts vector of values at the supporting points
     */
    PiecewiseConstantInterpolation(std::vector<double> supporting_points,
                                   std::vector<double>
                                       values_at_supp_pnts);

    /// \copydoc PiecewiseConstantInterpolation(std::vector<double>,
    /// std::vector<double>)
    PiecewiseConstantInterpolation(std::vector<int> const& supporting_points,
                                   std::vector<double>
                                       values_at_supp_pnts);

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
    double value(double pnt_to_interpolate) const;

private:
    std::vector<double> supp_pnts_;
    std::vector<double> values_at_supp_pnts_;
};
}  // end namespace MathLib
