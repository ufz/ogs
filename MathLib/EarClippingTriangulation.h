/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file EarClippingTriangulation.h
 *
 * Created on 2011-02-23 by Thomas Fischer
 */

#ifndef EARCLIPPINGTRIANGULATION_H_
#define EARCLIPPINGTRIANGULATION_H_

// STL
#include <list>

// GeoLib
#include "Polygon.h"
#include "Triangle.h"

// MathLib
#include "AnalyticalGeometry.h"

namespace MathLib {

class EarClippingTriangulation {
public:
	EarClippingTriangulation(const GeoLib::Polygon* ply, std::list<GeoLib::Triangle> &triangles, bool rot = true);
	virtual ~EarClippingTriangulation();
private:
	/**
	 * copies the points of the polygon to the vector _pnts
	 */
	inline void copyPolygonPoints (const GeoLib::Polygon* polygon);
	inline void rotate ();
	inline void ensureCWOrientation ();

	inline bool isEar(std::size_t v0, std::size_t v1, std::size_t v2) const;

	inline void initVertexList ();
	inline void initLists ();
	inline void clipEars ();

	/**
	 * a copy of the polygon points
	 */
	std::vector<GeoLib::Point*> _pnts;
	std::list<std::size_t> _vertex_list;
	std::list<std::size_t> _convex_vertex_list;
	std::list<std::size_t> _ear_list;

	/**
	 * triangles of the triangulation (maybe in the wrong orientation)
	 */
	std::list<GeoLib::Triangle> _triangles;

	MathLib::Orientation _original_orient;
};

} // end namespace MathLib

#endif /* EARCLIPPINGTRIANGULATION_H_ */
