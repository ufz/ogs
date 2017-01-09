/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   PiecewiseLinearMonotonicCurve.h
 *
 * Created on November 11, 2016, 10:49 AM
 */

#ifndef OGS_PIECE_WISE_LINEAR_MONOTONIC_CURVE_H
#define OGS_PIECE_WISE_LINEAR_MONOTONIC_CURVE_H

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MathLib
{
/// Class for strong monotonic curves
class PiecewiseLinearMonotonicCurve final : public PiecewiseLinearInterpolation
{
public:
    /**
     * @ x x coordinates of curve points
     * @ y y coordinates of curve points
     */
    PiecewiseLinearMonotonicCurve(std::vector<double>&& x,
                                  std::vector<double>&& y);

    ~PiecewiseLinearMonotonicCurve() = default;

    /// Get variable, x, or abscissa, by a given value \c y, the ordinate.
    /// If this curve is not monotonic, this function gives a fatal error.
    double getInversVariable(const double y) const;

private:
    /**
     *  Check monotonicity of this curve such as
     *  for any \f$i\f$,
     *    \f$y(x_{i+1})>y(x_{i})\f$ if \f$x_{i+1}>x_{i}\f$
     *  or for any \f$i\f$
     *    \f$y(x_{i+1})<y(x_{i})\f$ if \f$x_{i+1}>x_{i}\f$
     *  \return True if the curve exhibits strong monotonic.
     */
    bool isStrongMonotonic() const;
};
}  // namespace
#endif /* OGS_PIECE_WISE_LINEAR_MONOTONIC_CURVE_H */
