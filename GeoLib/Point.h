/***************************************************************************
 *  Point.h
 *  Created on: Jan 12, 2010
 *      Author: TF
**************************************************************************/

#ifndef POINT_H_
#define POINT_H_

#include "TemplatePoint.h"

namespace GEOLIB {

/**
 * \ingroup GEOLIB
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

}

/**
 * lexicographic comparison of points
 */
bool operator<= (GEOLIB::Point const & p0, GEOLIB::Point const & p1);

#endif /* POINT_H_ */
