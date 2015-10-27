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
#include "MathLib/LinAlg/Dense/DenseMatrix.h"

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
	MathLib::Vector3 n1(p0.getNormal());
	MathLib::Vector3 n2(p1.getNormal());
	MathLib::Vector3 direction(MathLib::crossProduct(n1, n2));

	// If the two planes p1 and p2 are in parallel the direction vector will be
	// zero. Either the planes have not any point in common or the planes are
	// identical.
	if (direction.getSqrLength() < eps*eps)
		return std::make_pair(direction, MathLib::Point3d());

	// determine a point p that fulfills both plane equations in Hessian normal
	// form (n1,p) - d1 = 0 and (n2,p) - d2 = 0

	// create right hand side of the system of linear equations
	std::array<double, 2> rhs{{p0.getDistance(), p1.getDistance()}};

	// copy vector components into the matrix of the system of linear equations
	MathLib::DenseMatrix<double, std::size_t> mat(2,3);
	mat(0,0) = n1[0];
	mat(0,1) = n1[1];
	mat(0,2) = n1[2];
	mat(1,0) = n2[0];
	mat(1,1) = n2[1];
	mat(1,2) = n2[2];

	// full pivot search
	std::size_t i(0), j(0);
	double temp_pivot(0.0);
	for (std::size_t r(0); r<mat.getNRows(); ++r) {
		for (std::size_t c(0); c<mat.getNCols(); ++c) {
			if (std::abs(mat(r,c)) > temp_pivot) {
				i = r;
				j = c;
				temp_pivot = std::abs(mat(i,j));
			}
		}
	}

	// exchange rows and cols to get the according to the absolute value biggest
	// entry to position (0,0)
	if (i != 0) { // if necessary exchange rows
		std::swap(mat(0,0), mat(i,0));
		std::swap(mat(0,1), mat(i,1));
		std::swap(mat(0,2), mat(i,2));
		std::swap(rhs[0], rhs[i]);
	}
	if (j != 0) { // if necessary exchange columns
		std::swap(mat(0,0), mat(0,j));
		std::swap(mat(1,0), mat(1,j));
		std::swap(mat(2,0), mat(2,j));
		std::swap(n1[0], n1[j]);
		std::swap(n2[0], n2[j]);
	}

	// eliminate the entry (1,0) and apply changes to the other entries of the
	// column and the right hand side
	long double const l(mat(1,0) / mat(0,0));
	mat(1,0) = 0.0;
	mat(1,1) -= l*mat(0,1);
	mat(1,2) -= l*mat(0,2);
	rhs[1] -= l*rhs[0];

	MathLib::Point3d p;
	if (j != 2 && std::abs(direction[2]) >= eps) {
		// direction vector is not in parallel with xy plane and it is save to
		// set p[2] zero
		p[2] = 0.0;
		// solve
		p[1] = rhs[1] / mat(1,1);
		p[0] = (rhs[0] - mat(0,1) * p[1])/ mat(0,0);
		// apply pivot exchanging to the components of the point
		std::swap(p[0], p[j]);
	} else if (j != 1 && std::abs(direction[1]) >= eps){
		// Direction vector is in parallel with xy plane and the y component is
		// not zero. Consequently, it is save to set p[1] to zero.
		p[1] = 0.0;
	} else {
		p[0] = 0.0;
	}

	return std::make_pair(direction, p);
}

}

#endif
