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

bool lessX (GeoLib::Point const & p0, GeoLib::Point const & p1)
{
	if (p0[0] <= p1[0]) return true;
	return false;
}

bool lessY (GeoLib::Point const & p0, GeoLib::Point const & p1)
{
	if (p0[1] <= p1[1]) return true;
	return false;
}

bool lessZ (GeoLib::Point const & p0, GeoLib::Point const & p1)
{
	if (p0[2] <= p1[2]) return true;
	return false;
}

bool operator<= (const GeoLib::Point& p0, const GeoLib::Point& p1)
{
	const double tol(std::numeric_limits<double>::epsilon());

	if (fabs (p0[0]-p1[0]) > tol * fabs(p0[0])) {
		if (p0[0] < p1[0]) return true;
		else return false;
	} else {
		// assume p0[0] == p1[0]
		if (fabs (p0[1]-p1[1]) > tol * fabs(p0[1])) {
			if (p0[1] < p1[1]) return true;
			else return false;
		} else {
			// assume p0[1] == p1[1] and p0[0] == p1[0]
			if (p0[2] < p1[2]) return true;
			else return false;
		}
	}
}

bool lessEq(const GeoLib::Point& p0, const GeoLib::Point& p1)
{
	const double tol(std::numeric_limits<double>::epsilon());

	// test a relative and an absolute criterion
	if (fabs(p0[0]-p1[0]) > tol * std::min(fabs(p1[0]), fabs(p0[0])) && fabs(p0[0]-p1[0]) > tol) {
		if (p0[0] <= p1[0]) return true;
		else return false;
	} else {
		// assume p0[0] == p1[0]
		if (fabs (p0[1]-p1[1]) > tol * fabs(p0[1]) && fabs(p0[1]-p1[1]) > tol) {
			if (p0[1] <= p1[1]) return true;
			else return false;
		} else {
			// assume p0[1] == p1[1] and p0[0] == p1[0]
			if (fabs (p0[2]-p1[2]) > tol * fabs(p0[2]) && fabs(p0[2]-p1[2]) > tol) {
				if (p0[2] <= p1[2]) return true;
				else return false;
			} else {
				return true;
			}
		}
	}
}


} // end namespace GeoLib
