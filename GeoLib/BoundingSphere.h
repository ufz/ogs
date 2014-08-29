/**
 * \file   Calculation of a minimum bounding sphere for a vector of points
 * \author Karsten Rink
 * \date   2014-07-11
 * \brief  Definition of the BoundingSphere class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#ifndef BOUNDINGSPHERE_H_
#define BOUNDINGSPHERE_H_

#include <vector>

#include "Vector3.h"
#include "Point.h"

namespace GeoLib
{

class BoundingSphere
{
public:
	/// Constructor using no points
	BoundingSphere();
	/// Copy constructor
	BoundingSphere(const BoundingSphere &sphere);
	/// Point-Sphere
	BoundingSphere(const GeoLib::Point &p);
	/// Constructor using center and radius
	BoundingSphere(const GeoLib::Point &p, double radius);
	/// Sphere through two points
	BoundingSphere(const GeoLib::Point &p, const GeoLib::Point &q);
	/// Sphere through three points
	BoundingSphere(const GeoLib::Point &p, const GeoLib::Point &q, const GeoLib::Point &);
	/// Sphere through four points
	BoundingSphere(const GeoLib::Point &p, const GeoLib::Point &q, const GeoLib::Point &R, const Point &S);
	/// Bounding sphere of n points
	BoundingSphere(const std::vector<GeoLib::Point*> &points);
	~BoundingSphere() {}

	/// Returns the center point of the sphere
	GeoLib::Point getCenter() const { return GeoLib::Point(_center.getCoords()); }

	/// Returns the radius of the sphere
	double getRadius() const {return _radius; }

	/// Returns the squared distance of a point from the sphere (for points within the sphere distance is negative)
	double sqrPointDist(const GeoLib::Point pnt) const;

	/// Creates n_points random points located on the surface of the sphere (useful for visualisation)
	std::vector<GeoLib::Point*>* getRandomSpherePoints(std::size_t n_points) const;

private:
	/**
	 * Recursive method for calculating a minimal bounding sphere for an arbitrary number of points.
	 * Algorithm based on Bernd Gärtner: Fast and Robust Smallest Enclosing Balls. ESA99, pages 325-338, 1999.
	 * Code based on "Smallest Enclosing Spheres" implementation by Nicolas Capens on flipcode's Developer Toolbox (www.flipcode.com)
	 */
	static BoundingSphere recurseCalculation(std::vector<GeoLib::Point*> sphere_points, std::size_t current_index, std::size_t n_points, std::size_t n_boundary_points);

	double _radius;
	MathLib::Vector3 _center;
};

} // namespace

#endif /* BOUNDINGSPHERE_H_ */
