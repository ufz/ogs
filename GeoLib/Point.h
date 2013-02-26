/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-12
 * \brief  Definition of the Point class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "TemplatePoint.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 */

template<typename T> class GeoPoint : public MathLib::TemplatePoint<T>, public GeoLib::GeoObject
{
public:
	GeoPoint(T x1, T x2, T x3) :
		MathLib::TemplatePoint<T>(x1, x2, x3), GeoLib::GeoObject()
	{}

	GeoPoint() :
		MathLib::TemplatePoint<T>(), GeoLib::GeoObject()
	{}

	GeoPoint (T const* x) :
		MathLib::TemplatePoint<T>(x), GeoObject()
	{}
};

typedef GeoLib::GeoPoint<double> Point;

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
