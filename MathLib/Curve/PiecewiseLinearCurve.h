/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   PiecewiseLinearCurve.h
 *
 * Created on November 11, 2016, 10:49 AM
 */

#ifndef OGS_PIECE_WISE_LINEAR_CURVE_H
#define OGS_PIECE_WISE_LINEAR_CURVE_H

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MathLib
{
class PiecewiseLinearCurve final : public PiecewiseLinearInterpolation
{
public:
    /**
     * @ x x coordinates of curve points
     * @ y y coordinates of curve points
     */
    PiecewiseLinearCurve(std::vector<double>&& x, std::vector<double>&& y)
        : PiecewiseLinearInterpolation(std::move(x), std::move(y), false)
    {
    }

    ~PiecewiseLinearCurve() = default;

    /**
     *  Check monotonicity of this curve such as
     *  for any \f$i\f$,
     *    \f$y(x_{i+1})>y(x_{i})\f$ if \f$x_{i+1}>x_{i}\f$
     *  or for any \f$i\f$
     *    \f$y(x_{i+1})<y(x_{i})\f$ if \f$x_{i+1}>x_{i}\f$
     *  \return True if the curve exhibits strong monotonic.
     */
    bool isStrongMonotonic() const;

    /// Get variable by a given value \c y.
    //  monotonicity must be check for the curve for using this function
    double getVariable(const double y) const;
};
}  // namespace
#endif /* OGS_PIECE_WISE_LINEAR_CURVE_H */
