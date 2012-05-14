/*
 * \file AnalyticalGeometry.h
 *
 *  Created on: Mar 17, 2010
 *      Author: TF
 */

#ifndef ANALYTICAL_GEOMETRY_H_
#define ANALYTICAL_GEOMETRY_H_

// MathLib
#include "Vector3.h"
// GeoLib
#include "Triangle.h"

namespace GeoLib {
	class Polyline;
}

namespace MathLib {

enum Orientation {
	CW = 1,
	CCW = 2,
	COLLINEAR = 3
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
Orientation getOrientation (const GeoLib::Point* p0, const GeoLib::Point* p1, const GeoLib::Point* p2);

/**
 * compute a supporting plane (represented by plane_normal and the value d) for the polygon
 * Let \f$n\f$ be the plane normal and \f$d\f$ a parameter. Then for all points \f$p \in R^3\f$ of the plane
 * it holds \f$ n \cdot p + d = 0\f$
 * @param pnts points of a closed polyline describing a polygon
 * @param plane_normal the normal of the plane the polygon is located in
 * @param d parameter from the plane equation
 */
void getNewellPlane (const std::vector<GeoLib::Point*>& pnts, MathLib::Vector &plane_normal, double& d);

/**
 *
 * @param plane_normal
 * @param pnts
 */
void rotatePointsToXY(MathLib::Vector &plane_normal, std::vector<GeoLib::Point*> &pnts);

bool isPointInTriangle (const GeoLib::Point* p,
		const GeoLib::Point* a, const GeoLib::Point* b, const GeoLib::Point* c);

/**
 * test for intersections of the line segments of the Polyline
 * @param ply the polyline
 * @param idx0 beginning index of the first line segment that has an intersection
 * @param idx1 beginning index of the second line segment that has an intersection
 * @param intersection_pnt the intersection point if the line segments intersect
 * @return true, if the polyline contains intersections
 */
bool lineSegmentsIntersect (const GeoLib::Polyline* ply, size_t &idx0, size_t &idx1, GeoLib::Point& intersection_pnt);

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

} // end namespace MathLib

#endif /* MATHTOOLS_H_ */
