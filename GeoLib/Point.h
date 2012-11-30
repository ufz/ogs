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
 * lexicographic comparison of points taking an epsilon into account
 */
bool lessEq(const GeoLib::Point& p0, const GeoLib::Point& p1);
}


#endif /* POINT_H_ */
