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

// GeoLib
#include "GeoObject.h"

// MathLib
#include "MathLib/Point3dWithID.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 */

// forward declaration
class PointVec;

class Point : public MathLib::Point3dWithID, public GeoLib::GeoObject
{
public:
	Point(double x1, double x2, double x3,
		std::size_t id = std::numeric_limits<std::size_t>::max()) :
		MathLib::Point3dWithID(std::array<double,3>({{x1, x2, x3}}), id),
		GeoLib::GeoObject()
	{}

	Point() :
		MathLib::Point3dWithID(), GeoLib::GeoObject()
	{}

	Point(double const* x) :
		MathLib::Point3dWithID(std::array<double,3>({{x[0], x[1], x[2]}}), 0),
		GeoLib::GeoObject()
	{}

	Point(Point const& p, std::size_t id) :
		MathLib::Point3dWithID(p, id), GeoLib::GeoObject()
	{}

	Point(std::array<double,3> const& x,
		std::size_t id = std::numeric_limits<std::size_t>::max()) :
		MathLib::Point3dWithID(x, id), GeoLib::GeoObject()
	{}

	/// return a geometry type
	virtual GEOTYPE getGeoType() const {return GEOTYPE::POINT;}

protected:
	friend GeoLib::PointVec;
	/// Resets the id.
	void setID(std::size_t id) { _id = id; }
};

static const Point ORIGIN(0, 0, 0);
}

#endif /* POINT_H_ */

