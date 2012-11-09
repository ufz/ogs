/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file PointWithID.h
 *
 *  Created on 2011-01-05 by Thomas Fischer
 */

#ifndef POINTWITHID_H_
#define POINTWITHID_H_

#include <limits>

#include "Point.h"

namespace GeoLib
{
/**
 * class PointWithID is derived from class Point in
 * order to extend the class Point with an ID.
 */
class PointWithID: public Point {
public:
	PointWithID(double x0, double x1, double x2, std::size_t id = std::numeric_limits<std::size_t>::max()) :
		Point(x0, x1, x2), _id(id)
	{}

	PointWithID(double const* const coords, std::size_t id = std::numeric_limits<std::size_t>::max()) :
		Point(coords), _id(id)
	{}

	PointWithID(GeoLib::Point const& pnt, std::size_t id) :
		Point(pnt), _id(id)
	{}

	/**
	 * standard constructor that initializes the id with 0 and calls
	 * the standard constructor of class Point
	 */
	PointWithID() :
		Point(), _id(0)
	{}

	std::size_t getID() const { return _id; }

protected:
	std::size_t _id;
};
} // end namespace GeoLib

#endif /* POINTWITHID_H_ */
