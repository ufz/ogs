/**
 * \file
 * \author Thomas Fischer
 * \date   2010-03-17
 * \brief  Definition of analytical geometry functions.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ANALYTICAL_GEOMETRY_H_
#define ANALYTICAL_GEOMETRY_H_

#include <memory>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/Vector3.h"

#include "Polygon.h"

namespace MathLib
{
    class PointVec;
}

namespace GeoLib
{
class Polyline;
class LineSegment;

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
 * it holds \f$ n \cdot p + d = 0\f$. The Newell algorithm is described in
 * \cite Ericson:2004:RCD:1121584 .
 * @param pnts_begin Iterator pointing to the initial point of a closed polyline describing a polygon
 * @param pnts_end Iterator pointing to the element following the last point of a closed polyline describing a polygon
 * @param plane_normal the normal of the plane the polygon is located in
 * @param d parameter from the plane equation
 */
template <typename InputIterator>
void getNewellPlane (InputIterator pnts_begin, InputIterator pnts_end,
                     MathLib::Vector3 &plane_normal,
                     double& d);

/**
 * compute a supporting plane (represented by plane_normal and the value d) for the polygon
 * Let \f$n\f$ be the plane normal and \f$d\f$ a parameter. Then for all points \f$p \in R^3\f$ of the plane
 * it holds \f$ n \cdot p + d = 0\f$. The Newell algorithm is described in
 * \cite Ericson:2004:RCD:1121584 .
 * @param pnts points of a closed polyline describing a polygon
 * @param plane_normal the normal of the plane the polygon is located in
 * @param d parameter from the plane equation
 */
template <class T_POINT>
void getNewellPlane (const std::vector<T_POINT*>& pnts,
                     MathLib::Vector3 &plane_normal,
                     double& d);

/** Same as getNewellPlane(pnts, plane_normal, d).
 */
template <class T_POINT>
std::pair<MathLib::Vector3, double> getNewellPlane(const std::vector<T_POINT>& pnts);

/**
 * Computes a rotation matrix that rotates the given 2D normal vector parallel to X-axis
 * @param v        a 2D normal vector to be rotated
 * @param rot_mat  a 2x2 rotation matrix
 */
template<class T_MATRIX>
void compute2DRotationMatrixToX(MathLib::Vector3 const& v, T_MATRIX & rot_mat);

/**
 * Computes a rotation matrix that rotates the given 3D normal vector parallel to X-axis
 * @param v        a 3D normal vector to be rotated
 * @param rot_mat  a 3x3 rotation matrix
 */
template <class T_MATRIX>
void compute3DRotationMatrixToX(MathLib::Vector3  const& v, T_MATRIX & rot_mat);

/**
 * Method computes the rotation matrix that rotates the given vector parallel to the \f$z\f$ axis.
 * @param plane_normal the (3d) vector that is rotated parallel to the \f$z\f$ axis
 * @param rot_mat 3x3 rotation matrix
 */
template <class T_MATRIX>
void computeRotationMatrixToXY(MathLib::Vector3 const& plane_normal,
                               T_MATRIX & rot_mat);

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
 * @param pnts_begin Iterator pointing to the initial element in a vector of points to be rotated
 * @param pnts_end Iterator pointing to the element following the last element in a vector of points to be rotated
 */
template <typename InputIterator>
void rotatePoints(
        MathLib::DenseMatrix<double> const& rot_mat,
        InputIterator pnts_begin, InputIterator pnts_end);

/**
 * rotate points according to the rotation matrix
 * @param rot_mat 3x3 dimensional rotation matrix
 * @param pnts vector of points
 */
template <typename P>
void rotatePoints(MathLib::DenseMatrix<double> const& rot_mat,
    std::vector<P*> const& pnts);

/**
 * rotate points to X-Y plane
 * @param pnts a vector of points with a minimum length of three.
 * Points are rotated using a rotation matrix computed from the first three points
 * in the vector. Point coordinates are modified as a result of the rotation.
 */
MathLib::DenseMatrix<double> rotatePointsToXY(std::vector<GeoLib::Point*> &pnts);

/**
 * rotate points to X-Y plane
 * @param p_pnts_begin Iterator pointing to the initial element in a vector of points used for computing a rotation matrix
 * @param p_pnts_end Iterator pointing to the element following the last point in a vector of points used for computing a rotation matrix
 * @param r_pnts_begin Iterator pointing to the initial element in a vector of points to be rotated
 * @param r_pnts_end Iterator pointing to the element following the last point in a vector of points to be rotated
 * Points are rotated using a rotation matrix computed from the first three points
 * in the vector. Point coordinates are modified as a result of the rotation.
 */
template <typename InputIterator1, typename InputIterator2>
MathLib::DenseMatrix<double> rotatePointsToXY(InputIterator1 p_pnts_begin,
                                              InputIterator1 p_pnts_end,
                                              InputIterator2 r_pnts_begin,
                                              InputIterator2 r_pnts_end);

/**
 * rotate points to X-Z plane
 * @param pnts a vector of points with a minimum length of three.
 * Points are rotated using a rotation matrix computed from the first three points
 * in the vector. Point coordinates are modified as a result of the rotation.
 */
void rotatePointsToXZ(std::vector<GeoLib::Point*> &pnts);

/**
 * test for intersections of the line segments of the Polyline
 * @param ply the polyline
 * @param seg_it0 iterator pointing to the first segment that has an intersection
 * @param seg_it1 iterator pointing to the second segment that has an intersection
 * @param intersection_pnt the intersection point if the segments intersect
 * @return true, if the polyline contains intersections
 */
bool lineSegmentsIntersect(const GeoLib::Polyline* ply,
                           GeoLib::Polyline::SegmentIterator& seg_it0,
                           GeoLib::Polyline::SegmentIterator& seg_it1,
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
 * @param s0 the first line segment.
 * @param s1 the second line segment.
 * @param s the intersection point if the segments do intersect
 * @return true, if the line segments intersect, else false
 */
bool lineSegmentIntersect(GeoLib::LineSegment const& s0,
                          GeoLib::LineSegment const& s1, GeoLib::Point& s);

/// A line segment is given by its two end-points. The function checks,
/// if the two line segments (ab) and (cd) intersects. This method checks the
/// intersection only in 2d.
/// @param ab first line segment
/// @param cd second line segment
/// @return empty vector in case there isn't an intersection point, a vector
/// containing one point if the line segments intersect in a single point, a
/// vector containing two points describing the line segment the original line
/// segments are interfering.
std::vector<MathLib::Point3d> lineSegmentIntersect2d(
    GeoLib::LineSegment const& ab, GeoLib::LineSegment const& cd);

/**
 * Calculates the intersection points of a line PQ and a triangle ABC.
 * This method requires ABC to be counterclockwise and PQ to point downward.
 * @return Intersection point or nullptr if there is no intersection.
 */
std::unique_ptr<GeoLib::Point> triangleLineIntersection(
    MathLib::Point3d const& a, MathLib::Point3d const& b,
    MathLib::Point3d const& c, MathLib::Point3d const& p,
    MathLib::Point3d const& q);

/**
 * Method first computes the intersection points of line segements of GeoLib::Polyline objects
 * (@see computeIntersectionPoints()) and pushes each intersection point in the GeoLib::PointVec
 * pnt_vec. For each intersection an id is returned.  This id is used to split the two
 * intersecting straight line segments in four straight line segments.
 */
void computeAndInsertAllIntersectionPoints(GeoLib::PointVec &pnt_vec,
    std::vector<GeoLib::Polyline*> & plys);

/**
 * Function rotates a polygon to the xy plane. For this reason, (1) the points of
 * the given polygon are copied, (2) a so called Newell plane is computed
 * (getNewellPlane()) and the points are rotated, (3) for security the
 * \f$z\f$ coordinates of the rotated points are set to zero and finally, (4) a
 * new polygon is constructed using the rotated points.
 * \see getNewellPlane()
 * @param polygon_in a copy of the polygon_in polygon will be rotated
 * @param plane_normal the normal of the original Newell plane
 * @return a rotated polygon
 */
GeoLib::Polygon rotatePolygonToXY(GeoLib::Polygon const& polygon_in,
    MathLib::Vector3 & plane_normal);

/// Sorts the vector of segments such that the \f$i\f$-th segment is connected
/// with the \f$i+1\f$st segment, i.e. the end point of the \f$i\f$-th segment
/// is the start point of the \f$i+1\f$st segment.
/// The current implementation requires that all segments have to be
/// connectable. In order to obtain a unique result the segments are sorted such
/// that the begin point of the first segment is \c seg_beg_pnt.
void sortSegments(MathLib::Point3d const& seg_beg_pnt,
    std::vector<GeoLib::LineSegment>& sub_segments);

} // end namespace GeoLib

#include "AnalyticalGeometry-impl.h"

#endif /* ANALYTICAL_GEOMETRY_H_ */
