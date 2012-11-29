/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file GeoMapper.h
 *
 * Created on 2012-09-25 by Karsten Rink
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

private:
	void mapData(MeshLib::Mesh const*const mesh = NULL);
	GeoLib::Grid<GeoLib::PointWithID>* getFlatGrid(MeshLib::Mesh const*const mesh, std::vector<GeoLib::PointWithID*> sfc_pnts) const;
	double getMeshElevation(double x, double y, MeshLib::Mesh const*const mesh) const;
	float getDemElevation(double x, double y) const;

	GeoLib::GEOObjects& _geo_objects;
	const std::string& _geo_name;

	// only necessary for mapping on mesh
	GeoLib::Grid<GeoLib::PointWithID>* _grid;

	// only necessary for mapping on DEM
	double _origin_x;
	double _origin_y;
	double _cellsize;
	unsigned _width;
	unsigned _height;
	float* _img_data;
};

#endif //GEOMAPPER_H
