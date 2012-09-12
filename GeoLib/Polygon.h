/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file Polygon.h
 *
 *  Created on 2010-06-21 by Thomas Fischer
 */

#ifndef POLYGON_H_
#define POLYGON_H_

// STL
#include <list>

// GeoLib
#include "AxisAlignedBoundingBox.h"
#include "Polyline.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 */

/**
 * edge classification
 */
class EdgeType
{
public:
	enum value {
		TOUCHING, //!< TOUCHING
		CROSSING, //!< CROSSING
		INESSENTIAL //!< INESSENTIAL
	};
};

/**
 *
 */
class Polygon : public Polyline
{
public:
	/**
	 * constructor checks if the given polyline is closed,
	 * and assures that the orientation is clock wise.
	 * @param ply closed Polyline
	 * @param init if true, check if polyline is closed, calculate bounding box
	 * @return
	 */
	Polygon(const Polyline &ply, bool init = true);

	Polygon (const std::vector<Point*>& pnt_vec);

	virtual ~Polygon();

	/**
	 *
	 * @return
	 */
	bool initialise ();

	/**
	 * Method checks if the given point is inside the polygon.
	 * The method requires that the polygon has clock wise orientation.
	 * @param pnt the Point
	 * @return if point is inside the polygon true, else false
	 */
	bool isPntInPolygon (const GeoLib::Point& pnt) const;
	/**
	 * wrapper for method isPntInPolygon (const GeoLib::Point&)
	 * @param x x coordinate of point
	 * @param y y coordinate of point
	 * @param z z coordinate of point
	 * @return if point is inside the polygon true, else false
	 */
	bool isPntInPolygon (double x, double y, double z) const;
	/**
	 * Method checks if all points of the polyline ply are inside of the polygon.
	 * @param ply the polyline that should be checked
	 * @return
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
	 * Calculates the next intersection point between the line segment (a,b) and the
	 * polygon starting with segment seg_num.
	 * @param a (input) the first point of the line segment
	 * @param b (input) the second point of the line segment
	 * @param intersection_pnt (output) next intersection point
	 * @param seg_num (input/output) the number of the polygon segment that is intersecting
	 */
	bool getNextIntersectionPointPolygonLine(GeoLib::Point const & a,
									GeoLib::Point const & b,
									GeoLib::Point* intersection_pnt,
									std::size_t& seg_num) const;


	void computeListOfSimplePolygons ();
	const std::list<Polygon*>& getListOfSimplePolygons ();

	friend bool operator==(Polygon const& lhs, Polygon const& rhs);
private:
	/**
	 * from book: Computational Geometry and Computer Graphics in C++, page 119
	 * get the type of edge with respect to the given point (2d method!)
	 * @param k number of line segment
	 * @param pnt point that is edge type computed for
	 * @return a value of enum EdgeType
	 */
	EdgeType::value getEdgeType (std::size_t k, GeoLib::Point const & pnt) const;

	void calculateAxisAlignedBoundingBox ();
	void ensureCWOrientation ();

	void splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it);
	void splitPolygonAtPoint (std::list<Polygon*>::iterator polygon_it);
	std::list<Polygon*> _simple_polygon_list;
	AABB _aabb;
};

/**
 * function creates a approximated circle area around a given point
 * @param middle_pnt the middle point of the circle
 * @param radius the radius of the circle
 * @param pnts (output) points that are used to approximate the circle
 * @param resolution number of point to use for approximation
 * @return a pointer to a polygon
 */
GeoLib::Polygon* createPolygonFromCircle (GeoLib::Point const& middle_pnt,
                                          double radius,
                                          std::vector<GeoLib::Point*> & pnts,
                                          std::size_t resolution = 12);

/**
 * comparison operator for polygons
 * @param lhs the first polygon
 * @param rhs the second polygon
 * @return true, if the polygons describe the same geometrical object
 */
bool operator==(Polygon const& lhs, Polygon const& rhs);

} // end namespace GeoLib

#endif /* POLYGON_H_ */
