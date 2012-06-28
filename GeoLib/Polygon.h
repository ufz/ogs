/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file Polygon.h
 *
 * Created on 2010-06-21 by Thomas Fischer
 */

#ifndef POLYGON_H_
#define POLYGON_H_

// STL
#include <list>

// GeoLib
#include "AxisAlignedBoundingBox.h"
#include "Polyline.h"

namespace GeoLib {

/**
 * \ingroup GeoLib
 */

/**
 * edge classification
 */
class EdgeType {
	public:
		enum value {
			TOUCHING,  //!< TOUCHING
			CROSSING,  //!< CROSSING
			INESSENTIAL//!< INESSENTIAL
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
	bool isPolylineInPolygon (const Polyline& ply) const;
	GeoLib::Point* getIntersectionPointPolygonLine (GeoLib::Point const & a, GeoLib::Point const & b) const;
	void computeListOfSimplePolygons ();
	const std::list<Polygon*>& getListOfSimplePolygons ();

private:
	/**
	 * get the type of edge with respect to the given point (2d method!)
	 * @param k number of line segment
	 * @param pnt point that is edge type computed for
	 * @return a value of enum EdgeType
	 */
	EdgeType::value getEdgeType (size_t k, GeoLib::Point const & pnt) const;

	void calculateAxisAlignedBoundingBox ();
	void ensureCWOrientation ();

	void splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it);
	void splitPolygonAtPoint (std::list<Polygon*>::iterator polygon_it);
	std::list<Polygon*> _simple_polygon_list;
	AABB _aabb;
};

GeoLib::Polygon* createPolygonFromCircle (GeoLib::Point const& middle_pnt, double radius,
		std::vector<GeoLib::Point*> & pnts, size_t resolution = 12);

} // end namespace GeoLib

#endif /* POLYGON_H_ */
