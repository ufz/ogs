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

namespace MeshGeoTools
{

MeshNodeSearcher::MeshNodeSearcher(MeshLib::Mesh const& mesh) :
		_mesh(mesh), _mesh_grid(_mesh.getNodes().cbegin(), _mesh.getNodes().cend())
{
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

	double epsilon_radius(0.001); // ToDo compute search length using the mesh
	// compute nodes (and supporting points) along polyline
	_mesh_nodes_along_polylines.push_back(
			MeshNodesAlongPolyline(_mesh.getNodes(), ply, epsilon_radius));
	return _mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size() - 1].getNodeIDs();
}


} // end namespace MeshGeoTools
