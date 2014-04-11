/**
 * \file   LayerVolumes.cpp
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Implementation of the LayerVolumes class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LayerVolumes.h"

#include <fstream>
#include <numeric>

#include "GEOObjects.h"
#include "PointVec.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "MeshGenerators/MeshLayerMapper.h"

LayerVolumes::LayerVolumes(GeoLib::GEOObjects &geo_objects)
: _invalid_value(-9999), _geo_objects(geo_objects)
{
}

bool LayerVolumes::createGeoVolumes(const MeshLib::Mesh &mesh, const std::vector<std::string> &raster_paths, double noDataReplacementValue)
{
	if (!allRastersExist(raster_paths))
		return false;

	std::vector<GeoLib::Point*> *points = new std::vector<GeoLib::Point*>;
	points->reserve(mesh.getNNodes() * raster_paths.size()); // max possible size
	std::vector<GeoLib::Surface*> *surfaces = new std::vector<GeoLib::Surface*>;
	surfaces->reserve(2 * raster_paths.size() - 1); // max possible size

	// DEM
	MeshLib::Mesh sfc_mesh (mesh);
	if (!MeshLayerMapper::LayerMapping(sfc_mesh, raster_paths[0], 0, 0, noDataReplacementValue))
		return false;

	const std::size_t nNodes (mesh.getNNodes());
	const std::vector<MeshLib::Node*> &nodes (mesh.getNodes());
	std::vector<std::size_t> pnts_above(nNodes);
	std::iota(pnts_above.begin(), pnts_above.end(), 0);
	std::vector<bool> node_status (nNodes, true);
	
	for (std::size_t i=0; i<nNodes; ++i)
		points->push_back(new GeoLib::Point(static_cast<GeoLib::Point>(*nodes[i])));
	surfaces->push_back(this->createSurface(mesh.getElements(), points, pnts_above, node_status));

	// subsurface layers
	const std::size_t nRasters (raster_paths.size());
	for (size_t i=1; i<nRasters; ++i)
	{
		if (!MeshLayerMapper::LayerMapping(sfc_mesh, raster_paths[0], 0, 0, _invalid_value))
		{
			cleanUpGeometryOnError(points, surfaces);
			return false;
		}
		this->addPoints(nodes, points, pnts_above, node_status);
		surfaces->push_back(this->createSurface(mesh.getElements(), points, pnts_above, node_status));
	}

	std::string geo_name (mesh.getName());
	_geo_objects.addPointVec(points, geo_name);
	const std::vector<std::size_t> id_map (_geo_objects.getPointVecObj(geo_name)->getIDMap());
	if (id_map.back() != id_map.size() - 1)
	{
		// TODO: fix indices for surfaces
	}
	_geo_objects.addSurfaceVec(surfaces, geo_name);

	// add boundary surfaces

	return true;
}

void LayerVolumes::addPoints(const std::vector<MeshLib::Node*> &nodes,
                             std::vector<GeoLib::Point*> *points,
                             std::vector<std::size_t> &pnts_above, 
                             std::vector<bool> &node_status) const
{
	std::size_t nNodes (nodes.size());
	for (std::size_t i=0; i<nNodes; ++i)
	{
		if ((*nodes[i])[2] == _invalid_value || (*nodes[i])[2] >= (*(*points)[pnts_above[i]])[2])
			node_status[i] = false;
		else
		{
			pnts_above[i] = points->size();
			node_status[i] = true;
			points->push_back(new GeoLib::Point(static_cast<GeoLib::Point>(*nodes[i])));
		}
	}
}

GeoLib::Surface* LayerVolumes::createSurface(const std::vector<MeshLib::Element*> &elements, 
                                             std::vector<GeoLib::Point*> *points,
                                             const std::vector<std::size_t> &pnts_above, 
                                             const std::vector<bool> &node_status) const
{
	const std::size_t nElems (elements.size());
	GeoLib::Surface* sfc = new GeoLib::Surface(*points);

	for (unsigned i=0; i<nElems; ++i)
	{
		MeshLib::Element* e (elements[i]);
		if (e->getGeomType() == MeshElemType::TRIANGLE &&
		   (node_status[e->getNodeIndex(0)] || node_status[e->getNodeIndex(1)] || node_status[e->getNodeIndex(2)]))
			sfc->addTriangle(pnts_above[e->getNodeIndex(0)], pnts_above[e->getNodeIndex(1)], pnts_above[e->getNodeIndex(2)]);
		if (e->getGeomType() == MeshElemType::QUAD)
		{
			if (node_status[e->getNodeIndex(0)] || node_status[e->getNodeIndex(1)] || node_status[e->getNodeIndex(2)])
				sfc->addTriangle(pnts_above[e->getNodeIndex(0)], pnts_above[e->getNodeIndex(1)], pnts_above[e->getNodeIndex(2)]);
			if (node_status[e->getNodeIndex(0)] || node_status[e->getNodeIndex(2)] || node_status[e->getNodeIndex(3)])
				sfc->addTriangle(pnts_above[e->getNodeIndex(0)], pnts_above[e->getNodeIndex(2)], pnts_above[e->getNodeIndex(3)]);
		}
		// all other element types are ignored (i.e. lines)
	}
	return sfc;
}

bool LayerVolumes::allRastersExist(const std::vector<std::string> &raster_paths) const
{
	for (auto raster = raster_paths.begin(); raster != raster_paths.end(); ++raster)
	{
		std::ifstream file_stream (*raster, std::ifstream::in);
		if (!file_stream.good())
			return false;
		file_stream.close();
	}
	return true;
}

void LayerVolumes::cleanUpGeometryOnError(std::vector<GeoLib::Point*> *points, std::vector<GeoLib::Surface*> *surfaces)
{
	std::for_each(points->begin(), points->end(), [](GeoLib::Point *point) { delete point; });
	std::for_each(surfaces->begin(), surfaces->end(), [](GeoLib::Surface *surface) { delete surface; });
	delete points;
	delete surfaces;
}
