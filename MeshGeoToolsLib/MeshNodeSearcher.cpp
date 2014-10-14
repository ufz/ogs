/**
 * @date Oct 24, 2013
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// MeshLib
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"

// MeshGeoToolsLib
#include "MeshGeoToolsLib/HeuristicSearchLength.h"
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"
#include "MeshGeoToolsLib/MeshNodesAlongSurface.h"

namespace MeshGeoToolsLib
{

MeshNodeSearcher::MeshNodeSearcher(MeshLib::Mesh const& mesh,
	MeshGeoToolsLib::SearchLength const& search_length_algorithm) :
		_mesh(mesh), _mesh_grid(_mesh.getNodes().cbegin(), _mesh.getNodes().cend()),
		_search_length(0.0)
{
	DBUG("Constructing MeshNodeSearcher obj.");
	//_search_length = search_length_algorithm.getSearchLength();

	double sum (0.0);
	double sum_of_sqr (0.0);
	const std::size_t ele_cnt(_mesh.getNElements());

	double min=0, max=0;
	for (const MeshLib::Element* e : _mesh.getElements()) {
		e->computeSqrNodeDistanceRange(min, max);
		sum += std::sqrt(min);
		sum_of_sqr += min;
	}

	const double mu (sum/ele_cnt);
	const double s (sqrt(1.0/(ele_cnt-1) * (sum_of_sqr - (sum*sum)/ele_cnt) ));
	// heuristic to prevent negative search lengths
	// in the case of a big standard deviation s
	double c(2.0);
	while (mu < c * s) {
		c *= 0.9;
	}

	_search_length = (mu - c * s)/2;

	DBUG("Calculated search length for mesh \"%s\" is %e.",
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
		boost::optional<std::size_t> node_id = this->getMeshNodeIDForPoint(*static_cast<const GeoLib::PointWithID*>(&geoObj));
		if (node_id) vec_nodes.push_back(*node_id);
		break;
	}
	case GeoLib::GEOTYPE::POLYLINE:
		vec_nodes = this->getMeshNodeIDsAlongPolyline(*static_cast<const GeoLib::Polyline*>(&geoObj));
		break;
	case GeoLib::GEOTYPE::SURFACE:
		vec_nodes = this->getMeshNodeIDsAlongSurface(*static_cast<const GeoLib::Surface*>(&geoObj));
		break;
	default:
		break;
	}
	return vec_nodes;
}

boost::optional<std::size_t> MeshNodeSearcher::getMeshNodeIDForPoint(GeoLib::Point const& pnt) const
{
	const MeshLib::Node* found = _mesh_grid.getNearestPoint(pnt.getCoords());
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
