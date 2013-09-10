 /**
 * \file
 * \author Karsten Rink
 * \date   2012-09-25
 * \brief  Definition of the GeoMapper class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

	void advancedMapOnMesh(const MeshLib::Mesh* mesh);

private:
	void mapData(MeshLib::Mesh const*const mesh = NULL);
	GeoLib::Grid<GeoLib::PointWithID>* getFlatGrid(MeshLib::Mesh const*const mesh, std::vector<GeoLib::PointWithID*> sfc_pnts) const;
	double getMeshElevation(double x, double y, MeshLib::Mesh const*const mesh) const;
	float getDemElevation(double x, double y) const;

	double getMaxSegmentLength(const std::vector<GeoLib::Polyline*> &lines) const;
	GeoLib::Point* calcIntersection(GeoLib::Point const*const p1, GeoLib::Point const*const p2, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const;
	bool isNodeOnLine(GeoLib::Point const*const p1, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const;
	bool isPntInBoundingBox(double ax, double ay, double bx, double by, double px, double py) const;
	std::size_t insertPointInLine(GeoLib::Polyline const*const line, GeoLib::Point const*const point, unsigned line_segment, const std::vector<unsigned> &line_segment_map) const;

	GeoLib::GEOObjects& _geo_objects;
	const std::string& _geo_name;

	// only necessary for mapping on mesh
	GeoLib::Grid<GeoLib::PointWithID>* _grid;

	// only necessary for mapping on DEM
	GeoLib::Raster *_raster;
};

#endif //GEOMAPPER_H
