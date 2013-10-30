/**
 * @file
 * @author git blame MeshNodeSearcher.cpp
 * @date Oct 24, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#include "Elements/Element.h"
#include "Elements/Line.h"

namespace MeshGeoTools
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
			it != elements.cend(); it++) {
		std::size_t const n_edges((*it)->getNEdges());
		for (std::size_t k(0); k<n_edges; k++) {
			MeshLib::Line const* line (dynamic_cast<MeshLib::Line const*>((*it)->getEdge(k)));
			double const len(line->getLength());
			sum += len;
			sum_of_sqr += len*len;
			delete line;
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
}

MeshNodeSearcher::~MeshNodeSearcher()
{
}

std::size_t MeshNodeSearcher::getMeshNodeIDForPoint(GeoLib::Point const& pnt) const
{
	return (_mesh_grid.getNearestPoint(pnt.getCoords()))->getID();
}

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsAlongPolyline(
		GeoLib::Polyline const& ply)
{
	std::vector<MeshNodesAlongPolyline>::const_iterator it(_mesh_nodes_along_polylines.begin());
	for (; it != _mesh_nodes_along_polylines.end(); it++) {
		if (it->getPolyline() == ply) {
			// we calculated mesh nodes for this polyline already
			return it->getNodeIDs();
		}
	}

	// compute nodes (and supporting points) along polyline
	_mesh_nodes_along_polylines.push_back(
			MeshNodesAlongPolyline(_mesh.getNodes(), ply, _search_length));
	return _mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size() - 1].getNodeIDs();
}

} // end namespace MeshGeoTools
