/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Point.h
 *
 * Created on 2010-01-12 by Thomas Fischer
 */

#ifndef POINT_H_
#define POINT_H_

#include "TemplatePoint.h"

namespace GeoLib {

/**
 * \ingroup GeoLib
 */

typedef TemplatePoint<double> Point;

/**
 * lexicographic comparison of points
 */
bool operator<= (GeoLib::Point const & p0, GeoLib::Point const & p1);

/**
 * lexicographical comparison of points taking an epsilon into account
 * @param p0 first input Point
 * @param p1 first input Point
 * @param tol tolerance (if in the comparison operation the property fabs(p0[k] - p1[k]) < tol
 * 	holds for the k-th coordinate the points are assumed the be equal in this coordinate)
 * @return true, if p0 is lexicographically smaller than p1
 */
bool lessEq(const GeoLib::Point& p0, const GeoLib::Point& p1, double tol = std::numeric_limits<double>::epsilon());
}


#endif /* POINT_H_ */
