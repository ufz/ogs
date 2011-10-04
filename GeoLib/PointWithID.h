/*
 * PointWithID.h
 *
 *  Created on: Jan 25, 2011
 *      Author: TF
 */

#ifndef POINTWITHID_H_
#define POINTWITHID_H_

#include "Point.h"

namespace GEOLIB {
/**
 * class PointWithID is derived from class Point in
 * order to extend the class Point with an ID.
 */
class PointWithID : public Point
{
public:
	PointWithID (double x0, double x1, double x2, size_t id) :
		Point (x0, x1, x2), _id (id)
	{}

	PointWithID (double const*const coords, size_t id) :
		Point (coords), _id (id)
	{}


	size_t getID () const { return _id; }

protected:
	size_t _id;
};
} // end namespace GEOLIB

#endif /* POINTWITHID_H_ */
