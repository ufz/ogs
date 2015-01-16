/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-12
 * \brief  Definition of the Point class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef POINT_H_
#define POINT_H_

// STL
#include <limits>

// GeoLib
#include "GeoObject.h"

// MathLib
#include "MathLib/MathPoint.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 */

class Point : public MathLib::MathPoint, public GeoLib::GeoObject
{
public:
	Point(double x1, double x2, double x3) :
		MathLib::MathPoint(), GeoLib::GeoObject()
	{
		this->_x[0] = x1;
		this->_x[1] = x2;
		this->_x[2] = x3;
	}

	Point() :
		MathLib::MathPoint(), GeoLib::GeoObject()
	{}

	Point (double const* x) :
		MathLib::MathPoint(x), GeoLib::GeoObject()
	{}

	/// return a geometry type
	virtual GEOTYPE getGeoType() const {return GEOTYPE::POINT;}
};

static const Point ORIGIN(0, 0, 0);

/**
 * lexicographic comparison of points
 */
bool operator<= (GeoLib::Point const & p0, GeoLib::Point const & p1);

/**
 * lexicographical comparison of points taking an epsilon into account
 * @param p0 first input Point
 * @param p1 first input Point
 * @param tol tolerance (if in the comparison operation the property fabs(p0[k] - p1[k]) < tol
 *     holds for the k-th coordinate the points are assumed the be equal in this coordinate)
 * @return true, if p0 is lexicographically smaller than p1
 */
bool lessEq(const GeoLib::Point& p0,
            const GeoLib::Point& p1,
            double tol = std::numeric_limits<double>::epsilon());
}

#endif /* POINT_H_ */
