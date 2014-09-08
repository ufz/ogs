 /**
 * \file
 * \author Karsten Rink
 * \date   2012-09-25
 * \brief  Definition of the GeoMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOMAPPER_H
#define GEOMAPPER_H

#include <cstddef>
#include <vector>

#include "GEOObjects.h"
#include "Point.h"
#include "Grid.h"

namespace MeshLib {
	class Mesh;
}

namespace GeoLib {
	class PointWithID;
	class Raster;
}

/**
 * \brief A set of tools for mapping the elevation of geometric objects
 */
class GeoMapper
{
public:
	GeoMapper(GeoLib::GEOObjects &geo_objects, const std::string &geo_name);
	~GeoMapper();

	void mapOnDEM(const std::string &file_name);
	void mapOnMesh(const std::string &file_name);
	void mapOnMesh(const MeshLib::Mesh* mesh);
	void advancedMapOnMesh(const MeshLib::Mesh* mesh, const std::string &new_geo_name);

private:
	// Manages the mapping geometric data (points, stations, boreholes) on a raster or mesh.
	void mapData();

	// Returns a grid containing all mesh surface points with elevation=0
	GeoLib::Grid<GeoLib::PointWithID>* getFlatGrid(MeshLib::Mesh const*const mesh, std::vector<GeoLib::PointWithID*> sfc_pnts) const;

	// Returns the elevation at Point (x,y) based on a mesh. This uses collision detection for triangles and nearest neighbor for quads.
	// NOTE: This medhod only returns correct values if the node numbering of the elements is correct!
	double getMeshElevation(double x, double y, double min_val, double max_val) const;

	// Returns the elevation at Point (x,y) based on a raster
	float getDemElevation(double x, double y) const;

	// Calculates the intersection of two lines embedded in the xy-plane
	GeoLib::Point* calcIntersection(GeoLib::Point const*const p1, GeoLib::Point const*const p2, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const;

	// Returns the position of a point within a line-segment
	unsigned getPointPosInLine(GeoLib::Polyline const*const line, unsigned start, unsigned end, GeoLib::Point const*const point, double eps) const;

	// Returns the maximum segment length in a polyline vector
	double getMaxSegmentLength(const std::vector<GeoLib::Polyline*> &lines) const;

	// Returns if a point p is within a bounding box defined by a and b
	bool isPntInBoundingBox(double ax, double ay, double bx, double by, double px, double py) const;

	GeoLib::GEOObjects& _geo_objects;
	std::string& _geo_name;

	// only necessary for mapping on mesh
	MeshLib::Mesh* _mesh;
	GeoLib::Grid<GeoLib::PointWithID>* _grid;

	// only necessary for mapping on DEM
	GeoLib::Raster *_raster;

};

#endif //GEOMAPPER_H
