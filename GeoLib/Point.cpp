/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file Point.cpp
 *
 * Created on 2010-06-22 by Thomas Fischer
 */


#include <cmath>
#include <limits>

#include "Point.h"

bool operator<= (const GeoLib::Point& p0, const GeoLib::Point& p1)
{
	double tol (sqrt (std::numeric_limits<double>::min()));

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


} // end namespace GeoLib
