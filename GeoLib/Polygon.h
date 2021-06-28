/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-21
 * \brief  Definition of the Polygon class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <list>
#include <vector>

#include "AABB.h"
#include "Point.h"
#include "Polyline.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 */

/**
 * edge classification
 */
enum class EdgeType
{
    TOUCHING, //!< TOUCHING
    CROSSING, //!< CROSSING
    INESSENTIAL //!< INESSENTIAL
};

/**
 * A polygon is a (closed) polyline. Thus class Polygon is derived from class Polyline.
 */
class Polygon : public Polyline
{
public:
    /**
     * constructor checks if the given polyline is closed,
     * and assures that the orientation is clock wise.
     * @param ply closed Polyline
     * @param init if true, check if polyline is closed, calculate bounding box
     */
    explicit Polygon(const Polyline& ply, bool init = true);

    Polygon(Polygon const& other);
    Polygon& operator=(Polygon const& rhs) = delete;

    ~Polygon() override;

    bool initialise ();

    /**
     * Method checks if the given point is inside the polygon.
     * The method requires that the polygon has clock wise orientation.
     * @param pnt the Point
     * @return if point is inside the polygon true, else false
     */
    bool isPntInPolygon(const MathLib::Point3d& pnt) const;

    /**
     * Checks if the straight line segment is contained within the polygon.
     * @param segment the straight line segment that is checked with
     * @return true if the straight line segment is within the polygon, else false
     */
    bool containsSegment(GeoLib::LineSegment const& segment) const;

    /**
     * Method checks if all points of the polyline ply are inside of the polygon.
     * @param ply the polyline that should be checked
     */
    bool isPolylineInPolygon (const Polyline& ply) const;
    /**
     * Method checks first if at least one (end!) point of a line segment of the polyline
     * is inside of the polygon. If this test fails each line segment of the polyline will
     * be tested against each polygon segment for intersection.
     * @param ply the polyline that should be checked
     * @return true if a part of the polyline is within the polygon
     */
    bool isPartOfPolylineInPolygon (const Polyline& ply) const;

    /**
     * Calculates the next intersection point between the line segment \c seg
     * and the polygon starting with segment \c seg_num.
     * @param seg (input) Line segment to compute intersection.
     * @param intersection_pnt (output) next intersection point
     * @param seg_num (in/out) the number of the polygon segment that is intersecting
     * @return true, if there was an intersection, i.e., the \c intersection_pnt
     * and \c seg_num contains new valid values
     */
    bool getNextIntersectionPointPolygonLine(GeoLib::LineSegment const& seg,
                                             GeoLib::Point& intersection_pnt,
                                             std::size_t& seg_num) const;

    friend bool operator==(Polygon const& lhs, Polygon const& rhs);
private:
    /**
     * Computes all intersections of the straight line segment and the polyline boundary
     * @param segment the line segment that will be processed
     * @return a possible empty vector containing the intersection points
     */
    std::vector<GeoLib::Point> getAllIntersectionPoints(
        GeoLib::LineSegment const& segment) const;

    /**
     * from book: Computational Geometry and Computer Graphics in C++, page 119
     * get the type of edge with respect to the given point (2d method!)
     * @param k number of line segment
     * @param pnt point that is edge type computed for
     * @return a value of enum EdgeType
     */
    EdgeType getEdgeType(std::size_t k, MathLib::Point3d const& pnt) const;

    void ensureCCWOrientation ();

    void splitPolygonAtIntersection(
        const std::list<Polygon*>::const_iterator& polygon_it);

    void splitPolygonAtPoint (const std::list<Polygon*>::iterator& polygon_it);
    std::list<Polygon*> _simple_polygon_list;
    AABB _aabb;
};

/**
 * comparison operator for polygons
 * @param lhs the first polygon
 * @param rhs the second polygon
 * @return true, if the polygons describe the same geometrical object
 */
bool operator==(Polygon const& lhs, Polygon const& rhs);

} // end namespace GeoLib
