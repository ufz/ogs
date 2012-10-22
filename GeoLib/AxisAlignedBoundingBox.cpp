/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file AxisAlignedBoundingBox.cpp
 *
 * Created on 2010-04-22 by Thomas Fischer
 */

#include "AxisAlignedBoundingBox.h"
#include <cmath>
#include <cstddef>
#include <limits>

namespace GeoLib
{
AABB::AABB()
{
	for (std::size_t k(0); k < 3; k++) {
		_min_pnt[k] = std::numeric_limits<double>::max();
		_max_pnt[k] = std::numeric_limits<double>::min();
	}
}

AABB::AABB(AABB const& src) :
	_min_pnt(src._min_pnt), _max_pnt(src._max_pnt)
{}

void AABB::update(GeoLib::Point const& pnt)
{
	for (size_t k(0); k<3; k++) {
		if (pnt[k] < _min_pnt[k])
			_min_pnt[k] = pnt[k];
		if (_max_pnt[k] < pnt[k])
			_max_pnt[k] = pnt[k];
	}
}

void AABB::update (double x, double y, double z)
{
	if (x < _min_pnt[0])
		_min_pnt[0] = x;
	if (_max_pnt[0] < x)
		_max_pnt[0] = x;
	if (y < _min_pnt[1])
		_min_pnt[1] = y;
	if (_max_pnt[1] < y)
		_max_pnt[1] = y;
	if (z < _min_pnt[2])
		_min_pnt[2] = z;
	if (_max_pnt[2] < z)
		_max_pnt[2] = z;
}

bool AABB::containsPoint (GeoLib::Point const & pnt, double eps) const
{
	return containsPoint (pnt[0], pnt[1], pnt[2], eps);
}

bool AABB::containsPoint (const double* pnt, double eps) const
{
	return containsPoint (pnt[0], pnt[1], pnt[2], eps);
}

bool AABB::containsPoint (double x, double y, double z, double eps) const
{
	if ((_min_pnt[0] <= x && x <= _max_pnt[0]) || std::fabs(_min_pnt[0] - x) < eps
					|| std::fabs(x - _max_pnt[0]) < eps) {
		if ((_min_pnt[1] <= y && y <= _max_pnt[1]) || std::fabs(_min_pnt[1] - y) < eps
						|| std::fabs(y - _max_pnt[1]) < eps) {
			if ((_min_pnt[2] <= z && z <= _max_pnt[2]) || std::fabs(_min_pnt[2] - z) < eps
							|| std::fabs(z - _max_pnt[2]) < eps) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	} else return false;
}

bool AABB::containsAABB (AABB const& other_aabb) const
{
	GeoLib::Point const& min_other(other_aabb.getMinPoint());
	GeoLib::Point const& max_other(other_aabb.getMaxPoint());
	for (unsigned k(0); k<3; k++) {
		if (_min_pnt[k] > min_other[k] || max_other[k] > _max_pnt[k])
			return false;
	}
	return true;
}

} // end namespace GeoLib
