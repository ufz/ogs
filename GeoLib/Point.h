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
 * comparison based on the x coordinate
 * @param p0 first point
 * @param p1 second point
 * @return true if the x coordinate of p0 is smaller equal the x coordinate of p1, else false
 */
bool lessX (Point const & p0, Point const & p1);

/**
 * comparison based on the y coordinate
 * @param p0 first point
 * @param p1 second point
 * @return true if the y coordinate of p0 is smaller equal the y coordinate of p1, else false
 */
bool lessY (Point const & p0, Point const & p1);

/**
 * comparison based on the z coordinate
 * @param p0 first point
 * @param p1 second point
 * @return true if the z coordinate of p0 is smaller equal the z coordinate of p1, else false
 */
bool lessZ (Point const & p0, Point const & p1);

/**
 * lexicographic comparison of points
 */
bool operator<= (GeoLib::Point const & p0, GeoLib::Point const & p1);
}


#endif /* POINT_H_ */
