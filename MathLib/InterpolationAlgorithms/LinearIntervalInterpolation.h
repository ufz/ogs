/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-22
 * \brief  Definition of the LinearIntervalInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/Error.h"

namespace MathLib {

/**
 * @brief Class (template) LinearIntervalInterpolation is a functional object performing
 * an interval mapping \f$f: [a,b] \to [c,d]\f$.
 *
 * Input numeric type has to be a floating point type and must behave well under the
 * operations addition, subtraction, multiplication and division. Let \f$a, b, c, d\f$
 * objects supporting the mentioned operations. Under the condition
 * \f$a \neq b\f$ an instance of the class computes a value within the interval
 * \f$[c, d]\f$, i.e., \f$f: [a,b] \to [c,d]\f$.
 */
template <typename NUMERIC_TYPE>
class LinearIntervalInterpolation {
public:
    /**
     * Constructor of class template for a linear map \f$y = m \cdot x + n\f$.
     * Under the prerequisite \f$a \neq b\f$ it initializes the coefficients
     * \f$m\f$ and \f$n\f$ in a correct way.
     * @param a first endpoint of the first interval
     * @param b second endpoint of the first interval
     * @param c first endpoint of the second interval
     * @param d second endpoint of the second interval
     */
    LinearIntervalInterpolation(NUMERIC_TYPE a, NUMERIC_TYPE b, NUMERIC_TYPE c, NUMERIC_TYPE d);
    /**
     * Method computes the value at point \f$x\f$ obtained by linear interpolation.
     * @param x the point the interpolation value is searched for
     * @return the interpolation value at point \f$x\f$
     */
    inline NUMERIC_TYPE operator() (NUMERIC_TYPE x) const;

private:
    /**
     * the slope of the linear map
     */
    NUMERIC_TYPE _m;
    /**
     * the offset of the linear map for \f$x\f$ equals zero
     */
    NUMERIC_TYPE _n;
};

template <typename NUMERIC_TYPE>
LinearIntervalInterpolation<NUMERIC_TYPE>::LinearIntervalInterpolation(NUMERIC_TYPE a, NUMERIC_TYPE b,
                NUMERIC_TYPE c, NUMERIC_TYPE d) :
    _m (d-c), _n(0.0)
{
    if (b == a) {
        OGS_FATAL("LinearIntervalInterpolation::LinearIntervalInterpolation: a == b, empty interval");
    }
    _m /= (b-a);
    _n = c - _m * a;
}

template <typename NUMERIC_TYPE>
inline NUMERIC_TYPE LinearIntervalInterpolation<NUMERIC_TYPE>::operator() (NUMERIC_TYPE x) const
{
    return _m * x + _n;
}

} // end namespace MathLib
