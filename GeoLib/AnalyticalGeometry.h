/**
 * \file
 * \author Thomas Fischer
 * \date   2010-03-17
 * \brief  Definition of analytical geometry functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ANALYTICAL_GEOMETRY_H_
#define ANALYTICAL_GEOMETRY_H_

// GeoLib
#include "Triangle.h"

// MathLib
#include "LinAlg/Dense/DenseMatrix.h"
#include "Vector3.h"

namespace GeoLib
{
class Polyline;

enum Orientation
{
	CW = 1, CCW = 2, COLLINEAR = 3
};

/**
 * computes the orientation of the three 2D-Points given by their coordinates
 * p0_x, p0_y, p1_x, p1_y, p2_x and p2_y
 * \returns CW (clockwise), CCW (counterclockwise) or COLLINEAR (points are on a line)
 */
Orientation getOrientation (const double& p0_x, const double& p0_y,
	const double& p1_x, const double& p1_y,
	const double& p2_x, const double& p2_y);

/**
 * wrapper for getOrientation ()
 */
Orientation getOrientation (const GeoLib::Point* p0,
                            const GeoLib::Point* p1,
                            const GeoLib::Point* p2);

/**
 * compute a supporting plane (represented by plane_normal and the value d) for the polygon
 * Let \f$n\f$ be the plane normal and \f$d\f$ a parameter. Then for all points \f$p \in R^3\f$ of the plane
 * it holds \f$ n \cdot p + d = 0\f$
 * @param pnts points of a closed polyline describing a polygon
 * @param plane_normal the normal of the plane the polygon is located in
 * @param d parameter from the plane equation
 */
void getNewellPlane (const std::vector<GeoLib::Point*>& pnts,
                     MathLib::Vector3 &plane_normal,
                     double& d);

/**
 * The vector plane_normal should be the surface normal of the plane surface described
 * by the points within the vector pnts. See function getNewellPlane() to get the
 * "plane normal" of a point set. The method rotates both the plane normal and
 * the points. The plane normal is rotated such that it is parallel to the \f$z\f$ axis.
 * @param plane_normal the normal of the plane
 * @param pnts pointers to points in a vector that should be rotated
 * @sa getNewellPlane()
 */
void rotatePointsToXY(MathLib::Vector3 &plane_normal, std::vector<GeoLib::Point*> &pnts);

/**
 * The vector plane_normal should be the surface normal of the plane surface described
 * by the points within the vector pnts. See function getNewellPlane() to get the
 * "plane normal" of a point set. The method rotates both the plane normal and
 * the points. The plane normal is rotated such that it is parallel to the \f$y\f$ axis.
 * @sa getNewellPlane()
 */
void rotatePointsToXZ(MathLib::Vector3 &plane_normal, std::vector<GeoLib::Point*> &pnts);

/**
 * Method computes the rotation matrix that rotates the given vector parallel to the \f$z\f$ axis.
 * @param plane_normal the (3d) vector that is rotated parallel to the \f$z\f$ axis
 * @param rot_mat 3x3 rotation matrix
 */
void computeRotationMatrixToXY(MathLib::Vector3 const& plane_normal,
                               MathLib::DenseMatrix<double> & rot_mat);

/**
 * rotate points according to the rotation matrix
 * @param rot_mat 3x3 dimensional rotation matrix
 * @param pnts vector of points
 */
void rotatePoints(MathLib::DenseMatrix<double> const& rot_mat, std::vector<GeoLib::Point*> &pnts);

bool isPointInTriangle (const GeoLib::Point* p,
		const GeoLib::Point* a, const GeoLib::Point* b, const GeoLib::Point* c);

bool isPointInTriangle(GeoLib::Point const& p,
				GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c,
				double eps = std::numeric_limits<double>::epsilon());

/**
 * test for intersections of the line segments of the Polyline
 * @param ply the polyline
 * @param idx0 beginning index of the first line segment that has an intersection
 * @param idx1 beginning index of the second line segment that has an intersection
 * @param intersection_pnt the intersection point if the line segments intersect
 * @return true, if the polyline contains intersections
 */
bool lineSegmentsIntersect (const GeoLib::Polyline* ply,
                            std::size_t &idx0,
                            std::size_t &idx1,
                            GeoLib::Point& intersection_pnt);

/**
 * A line segment is given by its two end-points. The function checks,
 * if the two line segments (ab) and (cd) intersects. Up to now only
 * 2D line segments are handled!
 * @param a first end-point of the first line segment
 * @param b second end-point of the first line segment
 * @param c first end-point of the second line segment
 * @param d second end-point of the second line segment
 * @param s the intersection point
 * @return true, if the line segments intersect, else false
 */
bool lineSegmentIntersect (const GeoLib::Point& a, const GeoLib::Point& b,
		const GeoLib::Point& c, const GeoLib::Point& d, GeoLib::Point& s);

/**
 * Calculates the intersection points of a line PQ and a triangle ABC. 
 * This method requires ABC to be counterclockwise and PQ to point downward.
 * @return Intersection point or NULL if there is no intersection.
 */
GeoLib::Point* triangleLineIntersection(GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c, GeoLib::Point const& p, GeoLib::Point const& q);

/// Calculates the scalar triple (u x v) . w
double scalarTriple(GeoLib::Point const& u, GeoLib::Point const& v, GeoLib::Point const& w);

/** 
 * Checks if a and b can be placed on a plane such that c and d lie on different sides of said plane.
 * (In 2D space this checks if c and d are on different sides of a line through a and b.)
 * @param a first point on plane
 * @param b second point on plane
 * @param c first point to test
 * @param d second point to test
 * @return true, if such a plane can be found, false otherwise 
 */
 bool dividedByPlane(const GeoLib::Point& a, const GeoLib::Point& b, 
	 const GeoLib::Point& c, const GeoLib::Point& d);

 /// Checks if the four given points are located on a plane.
 bool coplanar(const GeoLib::Point& a, const GeoLib::Point& b, 
	 const GeoLib::Point& c, const GeoLib::Point& d);


} // end namespace GeoLib

#endif /* ANALYTICAL_GEOMETRY_H_ */
