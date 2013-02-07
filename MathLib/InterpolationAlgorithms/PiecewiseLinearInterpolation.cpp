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
                                                           const std::vector<double>& values_at_supp_pnts) :
	_supp_pnts(supporting_points), _values_at_supp_pnts(values_at_supp_pnts)
{
	BaseLib::quicksort(_supp_pnts, 0, _supp_pnts.size(), _values_at_supp_pnts);
}

PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(const std::vector<double>& supporting_points,
                                                           const std::vector<double>& values_at_supp_pnts,
                                                           const std::vector<double>& points_to_interpolate,
                                                           std::vector<double>& values_at_interpol_pnts) :
	_supp_pnts(supporting_points), _values_at_supp_pnts(values_at_supp_pnts)
{
	values_at_interpol_pnts.clear();
	BaseLib::quicksort(_supp_pnts, 0, _supp_pnts.size(), _values_at_supp_pnts);
	for (std::size_t k(0); k < points_to_interpolate.size(); k++)
		values_at_interpol_pnts.push_back(this->getValue(points_to_interpolate[k]));
}

PiecewiseLinearInterpolation::~PiecewiseLinearInterpolation()
{}

double PiecewiseLinearInterpolation::getValue(double pnt_to_interpolate)
{
	// search interval that has the point inside
	std::size_t interval_idx(std::numeric_limits<std::size_t>::max());
	for (std::size_t k(1); k < _supp_pnts.size()
	     && interval_idx == std::numeric_limits<std::size_t>::max(); k++)
		if (_supp_pnts[k - 1] <= pnt_to_interpolate
		    && pnt_to_interpolate <= _supp_pnts[k])
			interval_idx = k - 1;

	if (interval_idx == std::numeric_limits<std::size_t>::max()) {
		const double dist_first(fabs(_supp_pnts[0] - pnt_to_interpolate));
		const double dist_last(fabs(_supp_pnts[_supp_pnts.size() - 1] - pnt_to_interpolate));
		if (dist_first < dist_last)
			interval_idx = 0;
		else
			interval_idx = _supp_pnts.size() - 2;
	}

	// compute linear interpolation polynom: y = m * x + n
	long double m((_values_at_supp_pnts[interval_idx + 1] - _values_at_supp_pnts[interval_idx])
	              / (_supp_pnts[interval_idx + 1] - _supp_pnts[interval_idx]));
	long double n(_values_at_supp_pnts[interval_idx] - m * _supp_pnts[interval_idx]);

	return m * pnt_to_interpolate + n;
}
} // end MathLib
