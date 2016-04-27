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

#include "PiecewiseLinearInterpolation.h"

#include "BaseLib/quicksort.h"

namespace MathLib
{
PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(
    std::vector<double>&& supporting_points,
    std::vector<double>&& values_at_supp_pnts,
    bool supp_pnts_sorted)
    : _supp_pnts(std::move(supporting_points)),
      _values_at_supp_pnts(std::move(values_at_supp_pnts))
{
    if (!supp_pnts_sorted) {
        BaseLib::quicksort<double, double>(_supp_pnts, static_cast<std::size_t> (0),
                                           _supp_pnts.size(), _values_at_supp_pnts);
    }
}

double PiecewiseLinearInterpolation::getValue(double pnt_to_interpolate) const
{
	// search interval that has the point inside
	std::size_t interval_idx(std::numeric_limits<std::size_t>::max());
	if (pnt_to_interpolate <= _supp_pnts.front()) {
		interval_idx = 0;
	} else {
		if (_supp_pnts.back() <= pnt_to_interpolate) {
			interval_idx = _supp_pnts.size() - 2;
		} else {
			auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(), pnt_to_interpolate));
			interval_idx = std::distance(_supp_pnts.begin(), it) - 1;
		}
	}
	if (pnt_to_interpolate > _supp_pnts.back()){
		return _values_at_supp_pnts[_supp_pnts.size() - 1]; 
	}
	if (pnt_to_interpolate < _supp_pnts.front()){
		return _values_at_supp_pnts[0];
	}
	// compute linear interpolation polynom: y = m * (x - support[i]) + value[i]
	const double m((_values_at_supp_pnts[interval_idx + 1] - _values_at_supp_pnts[interval_idx])
					/ (_supp_pnts[interval_idx + 1] - _supp_pnts[interval_idx]));

    return m * (pnt_to_interpolate - _supp_pnts[interval_idx]) + _values_at_supp_pnts[interval_idx];
}

double PiecewiseLinearInterpolation::getSlope(double pnt_to_interpolate) const
{
	std::size_t interval_idx(std::numeric_limits<std::size_t>::max());
	//double eps = 1e-9;
	if (pnt_to_interpolate <= _supp_pnts.front()) {
		interval_idx = 0;
	}
	else {
		if (_supp_pnts.back() <= pnt_to_interpolate) {
			interval_idx = _supp_pnts.size() - 2;
		}
		else {
			auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(), pnt_to_interpolate));
			interval_idx = std::distance(_supp_pnts.begin(), it) - 1;
		}
	}
	const double m((_values_at_supp_pnts[interval_idx + 1] - _values_at_supp_pnts[interval_idx])
		/ (_supp_pnts[interval_idx + 1] - _supp_pnts[interval_idx]));
	//double m_test = (getValue(pnt_to_interpolate + eps) - getValue(pnt_to_interpolate - eps)) / 2 / eps;
	return m;
}


double PiecewiseLinearInterpolation::PressureSaturationDependency(double pnt_to_interpolate, bool invert) const
{
	double dpds, dsdp, lim;
	dpds = GetCurveDerivative(pnt_to_interpolate);
	lim = (-1.0) / DBL_EPSILON;
	//
	if (invert)
	{ // Return dSw/dPc
		dsdp = 1.0 / dpds;
		if (dsdp < lim)
			dsdp = lim;
		return dsdp;
	}
	else
	{ // Return dPc/dSw
		if (dpds < lim)
			dpds = lim;
		return dpds;
	}

}

double PiecewiseLinearInterpolation::GetCurveDerivative(double pnt_to_interpolate) const
{
	std::size_t interval_max(std::numeric_limits<std::size_t>::max());//anz
	std::size_t interval_idx(std::numeric_limits<std::size_t>::max());//anz
	interval_max = _supp_pnts.size();
	if (pnt_to_interpolate <= _supp_pnts.front()) {
		interval_idx = 0;
	}
	else {
		if (_supp_pnts.back() <= pnt_to_interpolate) {
			interval_idx = _supp_pnts.size() - 2;
		}
		else {
			auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(), pnt_to_interpolate));
			interval_idx = std::distance(_supp_pnts.begin(), it) - 1;
		}
	}
	double sw = getValue(pnt_to_interpolate);
	if (sw < _values_at_supp_pnts[interval_max - 1]){
		sw = _values_at_supp_pnts[interval_max - 1];
		interval_idx = interval_max - 1;
	}
	else if (sw > _values_at_supp_pnts[0]){
		sw = _values_at_supp_pnts[0];
		interval_idx = 1;
	}
	//interval_idx = interval_max - 1 - interval_idx;
	if (interval_idx > 1 && interval_idx < interval_max - 2){
		double s1 = (0.5*_supp_pnts[interval_idx] - 0.5*_supp_pnts[interval_idx + 2]) / (0.5*_values_at_supp_pnts[interval_idx] - 0.5*_values_at_supp_pnts[interval_idx + 2]);
		double s2 = (0.5*_supp_pnts[interval_idx-1] - 0.5*_supp_pnts[interval_idx + 1]) / (0.5*_values_at_supp_pnts[interval_idx-1] - 0.5*_values_at_supp_pnts[interval_idx + 1]);
		double w = (sw - _values_at_supp_pnts[interval_idx + 1]) / (_values_at_supp_pnts[interval_idx] - _values_at_supp_pnts[interval_idx + 1]);
		return (1. - w) * s1 + w * s2;
	}
	else{
		if (std::fabs(_values_at_supp_pnts[interval_idx] - _values_at_supp_pnts[interval_idx + 1])>DBL_MIN)
			return (_supp_pnts[interval_idx] - _supp_pnts[interval_idx + 1]) / (_values_at_supp_pnts[interval_idx] - _values_at_supp_pnts[interval_idx + 1]);
		else
			return 1 / DBL_EPSILON;
	}
	
}

} // end MathLib
