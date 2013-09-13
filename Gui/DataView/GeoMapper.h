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
	void mapOnMesh(MeshLib::Mesh* mesh);

private:
	// Manages the mapping geometric data (points, stations, boreholes) on a raster or mesh.
	void mapData();
	// Returns a grid containing all mesh surface points with elevation=0
	GeoLib::Grid<GeoLib::PointWithID>* getFlatGrid(MeshLib::Mesh const*const mesh, std::vector<GeoLib::PointWithID*> sfc_pnts) const;
	// Returns the elevation at Point (x,y) based on a mesh. This uses collision detection for triangles and nearest neighbor for quads.
	double getMeshElevation(double x, double y, double min_val, double max_val) const;
	// Returns the elevation at Point (x,y) based on a raster
	float getDemElevation(double x, double y) const;

	GeoLib::GEOObjects& _geo_objects;
	const std::string& _geo_name;

	// only necessary for mapping on mesh
	MeshLib::Mesh* _mesh;
	GeoLib::Grid<GeoLib::PointWithID>* _grid;

	// only necessary for mapping on DEM
	GeoLib::Raster *_raster;

};

#endif //GEOMAPPER_H
