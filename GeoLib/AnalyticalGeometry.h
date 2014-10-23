/**
 * \file
 * \author Thomas Fischer
 * \date   2010-03-17
 * \brief  Definition of analytical geometry functions.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ANALYTICAL_GEOMETRY_H_
#define ANALYTICAL_GEOMETRY_H_

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/Vector3.h"

#include "Triangle.h"
#include "PointVec.h"


namespace GeoLib
{
class Polyline;

enum TriangleTest
{
	GAUSS, BARYCENTRIC
};

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
 * Method computes the rotation matrix that rotates the given vector parallel to the \f$z\f$ axis.
 * @param plane_normal the (3d) vector that is rotated parallel to the \f$z\f$ axis
 * @param rot_mat 3x3 rotation matrix
 */
void computeRotationMatrixToXY(MathLib::Vector3 const& plane_normal,
                               MathLib::DenseMatrix<double> & rot_mat);

/**
 * Method computes the rotation matrix that rotates the given vector parallel to the \f$y\f$ axis.
 * @param plane_normal the (3d) vector that is rotated parallel to the \f$y\f$ axis
 * @param rot_mat 3x3 rotation matrix
 */
void computeRotationMatrixToXZ(MathLib::Vector3 const& plane_normal,
                               MathLib::DenseMatrix<double> & rot_mat);

/**
 * rotate points according to the rotation matrix
 * @param rot_mat 3x3 dimensional rotation matrix
 * @param pnts vector of points
 */
void rotatePoints(MathLib::DenseMatrix<double> const& rot_mat, std::vector<GeoLib::Point*> &pnts);

/**
 * rotate points to X-Y plane
 * @param pnts a vector of points with a minimum length of three.
 * Points are rotated using a rotation matrix computed from the first three points
 * in the vector. Point coordinates are modified as a result of the rotation.
 */
void rotatePointsToXY(std::vector<GeoLib::Point*> &pnts);

/**
 * rotate points to X-Z plane
 * @param pnts a vector of points with a minimum length of three.
 * Points are rotated using a rotation matrix computed from the first three points
 * in the vector. Point coordinates are modified as a result of the rotation.
 */
void rotatePointsToXZ(std::vector<GeoLib::Point*> &pnts);

/**
 * Calculates the area of the triangle defined by its edge nodes a, b and c..
 * The formula is \f$A= \frac{1}{2} \cdot |u \times v|\f$, i.e. half of the area of the
 * parallelogram specified by the vectors\f$u=b-a\f$ and \f$v=c-a\f$.
 */
double calcTriangleArea(GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c);

/**
 * Calculates the volume of a tetrahedron.
 * The formula is V=1/6*|a(b x c)| with a=x1->x2, b=x1->x3 and c=x1->x4.
 */
double calcTetrahedronVolume(const double* x1, const double* x2, const double* x3, const double* x4);

/**
 * Tests if the given point p is within the triangle, defined by its edge nodes a, b and c.
 * Using two eps-values it is possible to test an 'epsilon' neighbourhood around the triangle
 * as well as an 'epsilon' outside the triangles plane.
 * @param p test point
 * @param a edge node of triangle
 * @param b edge node of triangle
 * @param c edge node of triangle
 * @param eps_pnt_out_of_plane eps allowing for p to be slightly off the plane spanned by abc
 * @param eps_pnt_out_of_tri eps allowing for p to be slightly off outside of abc
 * @param algorithm defines the method to use
 * @return true if the test point p is within the 'epsilon'-neighbourhood of the triangle
 */
bool isPointInTriangle(GeoLib::Point const& p,
                       GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c,
                       double eps_pnt_out_of_plane = std::numeric_limits<float>::epsilon(),
                       double eps_pnt_out_of_tri = std::numeric_limits<float>::epsilon(),
                       GeoLib::TriangleTest algorithm = GeoLib::GAUSS);

/**
 * Tests if the given point p is within the triangle, defined by its edge nodes a, b and c.
 * Using two eps-values it is possible to test an 'epsilon' neighbourhood around the triangle
 * as well as an 'epsilon' outside the triangles plane.
 * @param p test point
 * @param a edge node of triangle
 * @param b edge node of triangle
 * @param c edge node of triangle
 * @param eps_pnt_out_of_plane eps allowing for p to be slightly off the plane spanned by abc
 *                             ((orthogonal distance to the plane spaned by triangle)
 * @param eps_pnt_out_of_tri eps allowing for p to be slightly off outside of abc
 * @return true if the test point p is within the 'epsilon'-neighbourhood of the triangle
 */
bool gaussPointInTriangle(GeoLib::Point const& p,
                          GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c,
                          double eps_pnt_out_of_plane = std::numeric_limits<float>::epsilon(),
                          double eps_pnt_out_of_tri = std::numeric_limits<float>::epsilon());

/**
 * Tests if the given point p is within the triangle, defined by its edge nodes a, b and c.
 * Using two eps-values it is possible to test an 'epsilon' neighbourhood around the triangle
 * as well as an 'epsilon' outside the triangles plane.
 * Algorithm based on "Fundamentals of Computer Graphics" by Peter Shirley.
 * @param p test point
 * @param a edge node of triangle
 * @param b edge node of triangle
 * @param c edge node of triangle
 * @param eps_pnt_out_of_plane eps allowing for p to be slightly off the plane spanned by abc
 * @param eps_pnt_out_of_tri eps allowing for p to be slightly off outside of abc
 * @return true if the test point p is within the 'epsilon'-neighbourhood of the triangle
 */
bool barycentricPointInTriangle(GeoLib::Point const& p,
                                GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c,
                                double eps_pnt_out_of_plane = std::numeric_limits<float>::epsilon(),
                                double eps_pnt_out_of_tri = std::numeric_limits<float>::epsilon());

/**
 * Tests if the given point p is located within a tetrahedron spanned by points a, b, c, d.
 * If the tet specified by a, b, c, d is degenerated (i.e. all points are coplanar) the function
 * will return false because there is no meaningful concept of "inside" for such elements.
 * @param p test point
 * @param a edge node of tetrahedron
 * @param b edge node of tetrahedron
 * @param c edge node of tetrahedron
 * @param d edge node of tetrahedron
 * @param eps Accounts for numerical inaccuracies by allowing a point to be slightly outside of the element and still be regarded as inside.
 * @return true if the test point p is not located outside of abcd (i.e. inside or on a plane/edge).
 */
bool isPointInTetrahedron(GeoLib::Point const& p, GeoLib::Point const& a, GeoLib::Point const& b, 
                          GeoLib::Point const& c, GeoLib::Point const& d, double eps = std::numeric_limits<double>::epsilon());

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
 * Check if the two vectors \f$v, w \in R^3\f$ are in parallel
 * @param v first vector
 * @param w second vector
 * @return true if the vectors are in parallel, else false
*/
bool parallel(MathLib::Vector3 v, MathLib::Vector3 w);

/**
 * A line segment is given by its two end-points. The function checks,
 * if the two line segments (ab) and (cd) intersects.
 * @param a first end-point of the first line segment
 * @param b second end-point of the first line segment
 * @param c first end-point of the second line segment
 * @param d second end-point of the second line segment
 * @param s the intersection point if the segments do intersect
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
double scalarTriple(MathLib::Vector3 const& u, MathLib::Vector3 const& v, MathLib::Vector3 const& w);

/**
 * Checks if a point p is on the left or right side of a plane spanned by three points a, b, c.
 * @param p point to test
 * @param a first point on plane
 * @param b second point on plane
 * @param c third point on plane
 * @return If the triangle abc is ordered counterclockwise when viewed from p, the method will return a negative value,
 * otherwise it will return a positive value. If p is coplanar with abc, the function will return 0.
 */
double orientation3d(GeoLib::Point const& p,
                     GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c);

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
 bool isCoplanar(const GeoLib::Point& a, const GeoLib::Point& b, 
	 const GeoLib::Point& c, const GeoLib::Point& d);
 
/**
 * Method first computes the intersection points of line segements of GeoLib::Polyline objects
 * (@see computeIntersectionPoints()) and pushes each intersection point in the GeoLib::PointVec
 * pnt_vec. For each intersection an id is returned.  This id is used to split the two
 * intersecting straight line segments in four straight line segments.
 */
void computeAndInsertAllIntersectionPoints(
	GeoLib::PointVec &pnt_vec,
	std::vector<GeoLib::Polyline*> & plys);


} // end namespace GeoLib

#endif /* ANALYTICAL_GEOMETRY_H_ */
