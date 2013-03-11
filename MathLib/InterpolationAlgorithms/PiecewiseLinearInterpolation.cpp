/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-07
 * \brief  Implementation of the PiecewiseLinearInterpolation class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <cmath>
#include <limits>

#include "PiecewiseLinearInterpolation.h"

// BaseLib
#include "quicksort.h"

namespace MathLib
{
PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(const std::vector<double>& supporting_points,
                                                           const std::vector<double>& values_at_supp_pnts,
                                                           bool supp_pnts_sorted) :
	_supp_pnts(supporting_points), _values_at_supp_pnts(values_at_supp_pnts)
{
	if (!supp_pnts_sorted) {
		BaseLib::Quicksort<double, double>(_supp_pnts, static_cast<std::size_t> (0),
		                                   _supp_pnts.size(), _values_at_supp_pnts);
	}
}

PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(const std::vector<double>& supporting_points,
                                                           const std::vector<double>& values_at_supp_pnts,
                                                           const std::vector<double>& points_to_interpolate,
                                                           std::vector<double>& values_at_interpol_pnts,
                                                           bool supp_pnts_sorted) :
	_supp_pnts(supporting_points), _values_at_supp_pnts(values_at_supp_pnts)
{
	if (!supp_pnts_sorted) {
		BaseLib::Quicksort<double, double>(_supp_pnts, static_cast<std::size_t> (0),
		                                   _supp_pnts.size(),
		                                   _values_at_supp_pnts);
	}

	values_at_interpol_pnts.resize(points_to_interpolate.size());
	std::transform(points_to_interpolate.begin(), points_to_interpolate.end(),
	               values_at_interpol_pnts.begin(),
	               [&](double const& p) { return this->getValue(p); } );
}

PiecewiseLinearInterpolation::~PiecewiseLinearInterpolation()
{}

double PiecewiseLinearInterpolation::getValue(double pnt_to_interpolate) const
{
	// search interval that has the point inside
	std::size_t interval_idx(std::numeric_limits<std::size_t>::max());
	if (pnt_to_interpolate < _supp_pnts.front()) {
		interval_idx = 0;
	} else {
		if (_supp_pnts.back() < pnt_to_interpolate) {
			interval_idx = _supp_pnts.size() - 2;
		} else {
			auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(), pnt_to_interpolate));
			interval_idx = std::distance(_supp_pnts.begin(), it) - 1;
		}
	}

	// compute linear interpolation polynom: y = m * (x - support[i]) + value[i]
	const long double m((_values_at_supp_pnts[interval_idx + 1] - _values_at_supp_pnts[interval_idx])
					/ (_supp_pnts[interval_idx + 1] - _supp_pnts[interval_idx]));

	return m * (pnt_to_interpolate - _supp_pnts[interval_idx]) + _values_at_supp_pnts[interval_idx];
}
} // end MathLib
