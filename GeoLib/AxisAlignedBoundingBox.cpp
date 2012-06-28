/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
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
AABB::AABB ()
{
	for (std::size_t k(0); k < 3; k++)
	{
		_min_pnt[k] = std::numeric_limits<double>::max();
		_max_pnt[k] = std::numeric_limits<double>::min();
	}
}

AABB::AABB(AABB const& src) :
	_min_pnt(src._min_pnt.getCoords()), _max_pnt(src._max_pnt.getCoords())
{}

AABB::AABB ( const std::vector<GeoLib::Point*>* points )
{
	size_t nPoints (points->size());
	for (size_t i = 0; i < nPoints; i++)
		this->update((*(*points)[i])[0], (*(*points)[i])[1], (*(*points)[i])[2]);
}

void AABB::update (GeoLib::Point const & pnt)
{
	update (pnt[0], pnt[1], pnt[2]);
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

} // end namespace GeoLib
