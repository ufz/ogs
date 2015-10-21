/**
 * \brief  Definition of the class Plane - a plane in 3d space described in
 * Hessian normal form.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_PLANE_H
#define OGS_PLANE_H

#include <utility>

#include "MathLib/Vector3.h"

namespace GeoLib
{

class Plane final
{
public:
	/// Constructor that arguments already describe a plane in the hessian
	/// normal form.
	/// @param unit_normal normalized normal vector of the plane
	/// @attention{It is required that the plane normal is already normalized!}
	/// @param d the oriented distance to the parallel plane going through the
	/// origin
	Plane(MathLib::Vector3 const& unit_normal, double d)
		: _unit_normal(unit_normal), _oriented_distance(d)
	{}

	/// Constructs the hessian normal form of the plane spanned by the
	/// vectors \f$v\f$ and \f$w\f$ attached to the position vector \f$p\f$
	Plane(MathLib::Vector3 const& v, MathLib::Vector3 const& w,
		MathLib::Vector3 const& p)
		: _unit_normal(MathLib::crossProduct(v,w).getNormalizedVector()),
			_oriented_distance(-1.0 * MathLib::scalarProduct(_unit_normal, p))
	{}

	bool isPointInPlane(MathLib::Vector3 const& p,
		double tol=std::numeric_limits<float>::epsilon()) const
	{
		return getDistance(p)<tol;
	}

	/// Computes the distance between the position vector \f$p\f$ and the
	/// plane
	double getDistance(MathLib::Vector3 const& p) const
	{
		return MathLib::scalarProduct(_unit_normal,p)+_oriented_distance;
	}

	/// Returns the normalized normal of the plane.
	MathLib::Vector3 const& getNormal() const
	{
		return _unit_normal;
	}

	/// Returns the distance to the origin.
	double getDistance() const { return _oriented_distance; }

private:
	MathLib::Vector3 _unit_normal;
	double _oriented_distance;
};

/// Function computes the intersection of two planes.
/// @param p0 first plane
/// @param p1 second plane
/// @param eps tolerance value to check if vector components do not equal zero
/// @return a line described by the pair of a vector (direction) and a point on
/// the line
std::pair<MathLib::Vector3, MathLib::Point3d>
computePlanePlaneIntersection(GeoLib::Plane const& p0, GeoLib::Plane const& p1,
	double eps = std::numeric_limits<double>::epsilon())
{
	MathLib::Vector3 const& n1(p0.getNormal());
	MathLib::Vector3 const& n2(p1.getNormal());
	MathLib::Vector3 direction(MathLib::crossProduct(n1, n2));

	double d1(p0.getDistance());
	double d2(p1.getDistance());
	// determine a point p that fulfills both plane equations in Hessian normal
	// form (n1,p) - d1 = 0 and (n2,p) - d2 = 0
	MathLib::Point3d p;
	if (std::abs(direction[2]) >= eps) {
		// direction vector is not in parallel with xy plane and it is save to
		// set p[2] zero
		p[2] = 0.0;
		// solve the 2x3 system of linear equations manually
		double const l10(n2[0]/n1[0]);
		double const s(n2[1] - n1[1] * l10);
		p[1] = 1.0/s * (-d2 + d1*l10);
		p[0] = (-d1 - n1[1] * p[1])/n1[0];
	} else if (std::abs(direction[1]) >= eps){
		// Direction vector is in parallel with xy plane and the y component is
		// not zero. Consequently, it is save to set p[1] to zero.
		p[1] = 0.0;
		// solve the 2x3 system of linear equations manually
		double const l10(n2[0]/n1[0]);
		double const s(n2[2] - n1[2] * l10);
		p[2] = 1.0/s * (-d2 + d1*l10);
		p[0] = (-d1 - n1[2]*p[2])/n1[0];
	} else {
		// direction vector is parallel to e_x = (1,0,0)
		// => either (n1 = e_y and n2 = e_z) or (n1 = e_z and n2 = e_y)
		if (std::abs(n1[1]) >= eps) { // n1 = e_y and n2 = e_z
			p[1] = -d1 / n1[1];
			p[2] = -d2 / n2[2];
		} else { // n1 = e_z and n2 = e_y
			p[1] = -d2 / n2[1];
			p[2] = -d1 / n1[2];
		}
	}

	return std::make_pair(direction, p);
}

}

#endif
