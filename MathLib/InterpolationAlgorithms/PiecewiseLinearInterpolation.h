/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-07
 * \brief  Definition of the PiecewiseLinearInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PIECEWISELINEARINTERPOLATION_H_
#define PIECEWISELINEARINTERPOLATION_H_

#include <vector>

namespace MathLib
{
/**
 * This class implements a one dimensional piecewise linear interpolation algorithm.
 */
class PiecewiseLinearInterpolation final
{
public:
    /**
     * The constructor copies the entries of the vector of supporting points
     * \f$(x_0, x_1, \dots, x_n)\f$ and the entries of the vector of values at
     * the supporting points \f$(y_0, y_1, \dots, y_n)\f$ where \f$n\f$
     * is the number of entries of the vector. The number of supporting
     * points must be the same like the number of values at the supporting
     * points. It is assumed that \f$x_j\f$ corresponds to
     * \f$y_j\f$ for all \f$j \in [0, n]\f$.
     *
     * It is not assumed that the supporting points are sorted, i.e.
     * \f$x_0 < x_1 < \dots < x_n\f$. It is assumed, that the supporting points
     * are pairwise different. The user can set the flag supp_pnts_sorted to
     * true, if the supporting points are sorted. This will save some setup
     * time.
     * @param supporting_points vector of supporting points
     * @param values_at_supp_pnts vector of values at the supporting points
     * @param supp_pnts_sorted false (default), if it is sure the supporting points are sorted
     * one can set the switch to true
     */
    PiecewiseLinearInterpolation(std::vector<double>&& supporting_points,
                                 std::vector<double>&& values_at_supp_pnts,
                                 bool supp_pnts_sorted = false);

    /**
     * \brief Calculates the interpolation value.
     * @param pnt_to_interpolate The point should be located within the range
     * \f$[x_{\min}, x_{\max}]\f$, where \f$x_{\min} = \min_{1 \le j \le n} x_j\f$ and
     * \f$x_{\max} = \max_{1 \le j \le n} x_j\f$. Points outside of this interval are
     * extrapolated.
     * @return The interpolated value.
     */
    double getValue(double pnt_to_interpolate) const;

	double getSlope(double pnt_to_interpolate) const;

	double GetCurveDerivative(double pnt_to_interpolate) const;

	/*
	* return the derivative of the capillary pressure / saturation curve
	* if invert=true return dSw/dPc 
	* else if invert =false return dPc/dSw
	*/
	double PressureSaturationDependency(double pnt_to_interpolate, bool invert) const;
private:
    std::vector<double> _supp_pnts;
    std::vector<double> _values_at_supp_pnts;
};
} // end namespace MathLib

#endif /* PIECEWISELINEARINTERPOLATION_H_ */
