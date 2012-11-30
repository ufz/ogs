/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Point.cpp
 *
 * Created on 2010-06-22 by Thomas Fischer
 */


#include <cmath>
#include <limits>

#include "Point.h"

namespace GeoLib {

bool operator<= (const GeoLib::Point& p0, const GeoLib::Point& p1)
{
	if (p0[0] > p1[0]) {
		return false;
	} else {
		if (p0[0] < p1[0]) {
			return true;
		}
	}
	// => p0[0] == p1[0]

	if (p0[1] > p1[1]) {
		return false;
	} else {
		if (p0[1] < p1[1]) {
			return true;
		}
	}
	// => p0[1] == p1[1]

	if (p0[2] > p1[2]) {
		return false;
	} else {
		return true;
	}
}

bool lessEq(const GeoLib::Point& p0, const GeoLib::Point& p1, double tol)
{
	// test a relative and an absolute criterion
	if (fabs(p0[0]-p1[0]) > tol * std::min(fabs(p1[0]), fabs(p0[0])) && fabs(p0[0]-p1[0]) > tol) {
		if (p0[0] <= p1[0])
			return true;
		else
			return false;
	} else {
		// assume p0[0] == p1[0]
		if (fabs (p0[1]-p1[1]) > tol * fabs(p0[1]) && fabs(p0[1]-p1[1]) > tol) {
			if (p0[1] <= p1[1])
				return true;
			else
				return false;
		} else {
			// assume p0[1] == p1[1] and p0[0] == p1[0]
			if (fabs (p0[2]-p1[2]) > tol * fabs(p0[2]) && fabs(p0[2]-p1[2]) > tol) {
				if (p0[2] <= p1[2])
					return true;
				else
					return false;
			} else {
				return true;
			}
		}
	}
}


} // end namespace GeoLib
