/**
 * \file
 * \author Thomas Fischer
 * \date   2010-03-17
 * \brief  Definition of analytical geometry functions.
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "Polygon.h"

namespace GeoLib
{

enum Orientation
{
    CW = -1,
    COLLINEAR = 0,
    CCW = 1
};

/**
 * Computes the orientation of the three 2D-Points. This is a robust method.
 * \returns CW (clockwise), CCW (counterclockwise) or COLLINEAR (points are on a
 * line)
 */
Orientation getOrientation(MathLib::Point3d const& p0,
                           MathLib::Point3d const& p1,
                           MathLib::Point3d const& p2);
/**
 * Computes the orientation of the three 2D-Points. This is a non-robust method.
 * \returns CW (clockwise), CCW (counterclockwise) or COLLINEAR (points are on a
 * line)
 */
Orientation getOrientationFast(MathLib::Point3d const& p0,
                               MathLib::Point3d const& p1,
                               MathLib::Point3d const& p2);
/**
 * Compute a supporting plane (represented by plane_normal and the value d) for
 * the polygon.
 * Let \f$n\f$ be the plane normal and \f$d\f$ a parameter. Then for all points
 * \f$p \in R^3\f$ of the plane it holds \f$ n \cdot p + d = 0\f$.
 * The Newell algorithm is described in \cite Ericson:2004:RCD:1121584 .
 * \param pnts_begin Iterator pointing to the initial point of a closed polyline
 * describing a polygon
 * \param pnts_end Iterator pointing to the element following the last point of
 * a closed polyline describing a polygon
 * \return pair of plane_normal and the parameter d: the normal of the plane the
 * polygon is located in, d parameter from the plane equation
 */
template <typename InputIterator>
std::pair<Eigen::Vector3d, double> getNewellPlane(InputIterator pnts_begin,
                                                  InputIterator pnts_end);

/**
 * Compute a supporting plane (represented by plane_normal and the value d) for
 * the polygon.
 * Let \f$n\f$ be the plane normal and \f$d\f$ a parameter. Then for all points
 * \f$p \in R^3\f$ of the plane it holds \f$ n \cdot p + d = 0\f$.
 * The Newell algorithm is described in \cite Ericson:2004:RCD:1121584 .
 * \param pnts points of a closed polyline describing a polygon
 * \return pair of plane_normal and the parameter d: the normal of the plane the
 * polygon is located in, parameter d from the plane equation
 */
template <class T_POINT>
std::pair<Eigen::Vector3d, double> getNewellPlane(
    const std::vector<T_POINT*>& pnts);

/** Same as getNewellPlane(pnts).
 */
template <class T_POINT>
std::pair<Eigen::Vector3d, double> getNewellPlane(
    const std::vector<T_POINT>& pnts);

/**
 * Computes a rotation matrix that rotates the given 2D normal vector parallel
 * to X-axis
 * \param v a 2D normal vector to be rotated
 * \return a 3x3 rotation matrix where the upper, left, 2x2 block
 * contains the entries necessary for the 2D rotation
 */
Eigen::Matrix3d compute2DRotationMatrixToX(Eigen::Vector3d const& v);

/**
 * Computes a rotation matrix that rotates the given 3D normal vector parallel
 * to X-axis.
 * \param v        a 3D normal vector to be rotated
 * \return a 3x3 rotation matrix
 */
Eigen::Matrix3d compute3DRotationMatrixToX(Eigen::Vector3d const& v);

/**
 * Method computes the rotation matrix that rotates the given vector parallel to
 * the \f$z\f$ axis.
 * \param n the (3d) vector that is rotated parallel to the \f$z\f$ axis
 * \return rot_mat 3x3 rotation matrix
 */
Eigen::Matrix3d computeRotationMatrixToXY(Eigen::Vector3d const& n);

/**
 * rotate points according to the rotation matrix
 * \param rot_mat 3x3 dimensional rotation matrix
 * \param pnts_begin Iterator pointing to the initial element in a vector of
 * points to be rotated
 * \param pnts_end Iterator pointing to the element following the last element
 * in a vector of points to be rotated
 */
template <typename InputIterator>
void rotatePoints(Eigen::Matrix3d const& rot_mat, InputIterator pnts_begin,
                  InputIterator pnts_end);

/**
 * rotate points according to the rotation matrix
 * \param rot_mat 3x3 dimensional rotation matrix
 * \param pnts vector of points
 */
template <typename P>
void rotatePoints(Eigen::Matrix3d const& rot_mat, std::vector<P*> const& pnts);

/**
 * rotate points to X-Y plane
 * \param pnts a vector of points with a minimum length of three.
 * Points are rotated using a rotation matrix computed from the first three
 * points in the vector. Point coordinates are modified as a result of the
 * rotation.
 */
Eigen::Matrix3d rotatePointsToXY(std::vector<Point*>& pnts);

/**
 * rotate points to X-Y plane
 * \param p_pnts_begin Iterator pointing to the initial element in a vector of
 * points used for computing a rotation matrix
 * \param p_pnts_end Iterator pointing to the element following the last point
 * in a vector of points used for computing a rotation matrix
 * \param r_pnts_begin Iterator pointing to the initial element in a vector of
 * points to be rotated
 * \param r_pnts_end Iterator pointing to the element following the last point
 * in a vector of points to be rotated Points are rotated using a rotation
 * matrix computed from the first three points in the vector. Point coordinates
 * are modified as a result of the rotation.
 */
template <typename InputIterator1, typename InputIterator2>
Eigen::Matrix3d rotatePointsToXY(InputIterator1 p_pnts_begin,
                                 InputIterator1 p_pnts_end,
                                 InputIterator2 r_pnts_begin,
                                 InputIterator2 r_pnts_end);

/**
 * test for intersections of the line segments of the Polyline
 * \param ply the polyline
 * \param seg_it0 iterator pointing to the first segment that has an
 * intersection
 * \param seg_it1 iterator pointing to the second segment that has an
 * intersection
 * \param intersection_pnt the intersection point if the segments intersect
 * \return true, if the polyline contains intersections
 */
bool lineSegmentsIntersect(const Polyline* ply,
                           Polyline::SegmentIterator& seg_it0,
                           Polyline::SegmentIterator& seg_it1,
                           Point& intersection_pnt);

/**
 * Check if the two vectors \f$v, w \in R^3\f$ are in parallel
 * \param v first vector
 * \param w second vector
 * \return true if the vectors are in parallel, else false
 */
bool parallel(Eigen::Vector3d v, Eigen::Vector3d w);

/**
 * A line segment is given by its two end-points. The function checks,
 * if the two line segments (ab) and (cd) intersects.
 * \param s0 the first line segment.
 * \param s1 the second line segment.
 * \param s the intersection point if the segments do intersect
 * \return true, if the line segments intersect, else false
 */
bool lineSegmentIntersect(LineSegment const& s0, LineSegment const& s1,
                          Point& s);

/// A line segment is given by its two end-points. The function checks,
/// if the two line segments (ab) and (cd) intersects. This method checks the
/// intersection only in 2d.
/// @param ab first line segment
/// @param cd second line segment
/// @return empty vector in case there isn't an intersection point, a vector
/// containing one point if the line segments intersect in a single point, a
/// vector containing two points describing the line segment the original line
/// segments are interfering.
std::vector<MathLib::Point3d> lineSegmentIntersect2d(LineSegment const& ab,
                                                     LineSegment const& cd);

/**
 * Calculates the intersection points of a line PQ and a triangle ABC.
 * This method requires ABC to be counterclockwise and PQ to point downward.
 * \return Intersection point or nullptr if there is no intersection.
 */
std::unique_ptr<Point> triangleLineIntersection(MathLib::Point3d const& a,
                                                MathLib::Point3d const& b,
                                                MathLib::Point3d const& c,
                                                MathLib::Point3d const& p,
                                                MathLib::Point3d const& q);

/**
 * Method first computes the intersection points of line segments of Polyline
 * objects
 * (@see computeIntersectionPoints()) and pushes each intersection point in the
 * PointVec pnt_vec. For each intersection an id is returned.  This id is used
 * to split the two intersecting straight line segments in four straight line
 * segments.
 */
void computeAndInsertAllIntersectionPoints(PointVec& pnt_vec,
                                           std::vector<Polyline*>& plys);

/**
 * Function rotates a polygon to the xy plane. For this reason, (1) the points
 * of the given polygon are copied, (2) a so called Newell plane is computed
 * (getNewellPlane()) and the points are rotated, (3) for accuracy reasons the
 * \f$z\f$ coordinates of the rotated points are set to zero
 * \see getNewellPlane()
 * \param polygon_in a copy of the polygon_in polygon will be rotated
 * \return vector of rotated points and normal based on the original Newell
 * plane
 */
std::tuple<std::vector<GeoLib::Point*>, Eigen::Vector3d>
rotatePolygonPointsToXY(GeoLib::Polygon const& polygon_in);

/// Sorts the vector of segments such that the \f$i\f$-th segment is connected
/// with the \f$i+1\f$st segment, i.e. the end point of the \f$i\f$-th segment
/// is the start point of the \f$i+1\f$st segment.
/// The current implementation requires that all segments have to be
/// connectable. In order to obtain a unique result the segments are sorted such
/// that the begin point of the first segment is \c seg_beg_pnt.
void sortSegments(MathLib::Point3d const& seg_beg_pnt,
                  std::vector<LineSegment>& sub_segments);

}  // end namespace GeoLib

#include "AnalyticalGeometry-impl.h"
