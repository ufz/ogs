/**
 * @file
 * @date Oct 24, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

// MeshLib
#include "Elements/Element.h"
#include "Elements/Line.h"

// MeshGeoToolsLib
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"
#include "MeshGeoToolsLib/MeshNodesAlongSurface.h"

namespace MeshGeoToolsLib
{

MeshNodeSearcher::MeshNodeSearcher(MeshLib::Mesh const& mesh) :
		_mesh(mesh), _mesh_grid(_mesh.getNodes().cbegin(), _mesh.getNodes().cend()),
		_search_length(0.0)
{
	double sum (0.0);
	double sum_of_sqr (0.0);
	std::size_t edge_cnt(0);
	std::vector<MeshLib::Element*> const& elements(_mesh.getElements());

	for (std::vector<MeshLib::Element*>::const_iterator it(elements.cbegin());
			it != elements.cend(); ++it) {
		std::size_t const n_edges((*it)->getNEdges());
		for (std::size_t k(0); k<n_edges; k++) {
			MeshLib::Line const* edge(dynamic_cast<MeshLib::Line const*>((*it)->getEdge(k)));
			if (!edge) {
				delete edge;
				continue;
			}
			double const len(edge->getLength());
			sum += len;
			sum_of_sqr += len*len;
			delete edge;
		}
		edge_cnt += n_edges;
	}

	const double mu (sum/edge_cnt);
	const double s (sqrt(1.0/(edge_cnt-1) * (sum_of_sqr - (sum*sum)/edge_cnt) ));
	// heuristic to prevent negative search lengths
	// in the case of a big standard deviation s
	double c(2.0);
	while (mu < c * s) {
		c *= 0.9;
	}

	_search_length = (mu - c * s)/2;

	DBUG("[MeshNodeSearcher::MeshNodeSearcher] Calculated search length for mesh \"%s\" is %f.",
			_mesh.getName().c_str(), _search_length);
}

MeshNodeSearcher::~MeshNodeSearcher()
{
	std::vector<MeshNodesAlongPolyline*>::iterator it(_mesh_nodes_along_polylines.begin());
	for (; it != _mesh_nodes_along_polylines.end(); ++it) {
		delete (*it);
	}
}

std::vector<std::size_t> MeshNodeSearcher::getMeshNodeIDs(GeoLib::GeoObject const& geoObj)
{
	std::vector<std::size_t> vec_nodes;
	switch (geoObj.getGeoType()) {
	case GeoLib::GEOTYPE::POINT:
	{
		auto node_id = this->getMeshNodeIDForPoint(*dynamic_cast<const GeoLib::PointWithID*>(&geoObj));
		if (node_id) vec_nodes.push_back(*node_id);
		break;
	}
	case GeoLib::GEOTYPE::POLYLINE:
		vec_nodes = this->getMeshNodeIDsAlongPolyline(*dynamic_cast<const GeoLib::Polyline*>(&geoObj));
		break;
	case GeoLib::GEOTYPE::SURFACE:
		vec_nodes = this->getMeshNodeIDsAlongSurface(*dynamic_cast<const GeoLib::Surface*>(&geoObj));
		break;
	default:
		break;
	}
	return vec_nodes;
}

boost::optional<std::size_t> MeshNodeSearcher::getMeshNodeIDForPoint(GeoLib::Point const& pnt) const
{
	auto* found = _mesh_grid.getNearestPoint(pnt.getCoords());
	if (found)
		return found->getID();
	else
		return boost::none;
}

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsAlongPolyline(
		GeoLib::Polyline const& ply)
{
	return getMeshNodesAlongPolyline(ply).getNodeIDs();
}

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsAlongSurface(GeoLib::Surface const& sfc)
{
	return getMeshNodesAlongSurface(sfc).getNodeIDs();
}

MeshNodesAlongPolyline& MeshNodeSearcher::getMeshNodesAlongPolyline(GeoLib::Polyline const& ply)
{
	std::vector<MeshNodesAlongPolyline*>::const_iterator it(_mesh_nodes_along_polylines.begin());
	for (; it != _mesh_nodes_along_polylines.end(); ++it) {
		if (&(*it)->getPolyline() == &ply) {
			// we calculated mesh nodes for this polyline already
			return *(*it);
		}
	}

	// compute nodes (and supporting points) along polyline
	_mesh_nodes_along_polylines.push_back(
			new MeshNodesAlongPolyline(_mesh, ply, _search_length));
	return *_mesh_nodes_along_polylines.back();
}

MeshNodesAlongSurface& MeshNodeSearcher::getMeshNodesAlongSurface(GeoLib::Surface const& sfc)
{
	std::vector<MeshNodesAlongSurface*>::const_iterator it(_mesh_nodes_along_surfaces.begin());
	for (; it != _mesh_nodes_along_surfaces.end(); ++it) {
		if (&(*it)->getSurface() == &sfc) {
			// we calculated mesh nodes for this polyline already
			return *(*it);
		}
	}

	// compute nodes (and supporting points) along polyline
	_mesh_nodes_along_surfaces.push_back(
			new MeshNodesAlongSurface(_mesh, sfc));
	return *_mesh_nodes_along_surfaces.back();
}

} // end namespace MeshGeoTools
