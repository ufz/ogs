/**
 * @file
 * @date Aug 9, 2010
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// BaseLib
#include "quicksort.h"

// MeshGeoToolsLib
#include "MeshNodesAlongPolyline.h"

#include "MathTools.h"

#include <algorithm>

namespace MeshGeoToolsLib
{
MeshNodesAlongPolyline::MeshNodesAlongPolyline(
		MeshLib::Mesh const& mesh,
		GeoLib::Polyline const& ply,
		double epsilon_radius) :
	_mesh(mesh), _ply(ply)
{
	assert(epsilon_radius > 0);
	auto &mesh_nodes = _mesh.getNodes();
	const std::size_t n_nodes (mesh_nodes.size());
	// loop over all nodes
	for (size_t i = 0; i < n_nodes; i++) {
		double dist = _ply.getDistanceAlongPolyline(*mesh_nodes[i], epsilon_radius);
		if (dist >= 0.0) {
			_msh_node_ids.push_back(mesh_nodes[i]->getID());
			_dist_of_proj_node_from_ply_start.push_back(dist);
		}
	}

	// sort the nodes along the polyline according to their distances
	BaseLib::Quicksort<double> (_dist_of_proj_node_from_ply_start, 0,
					_dist_of_proj_node_from_ply_start.size(), _msh_node_ids);
}

MeshLib::Mesh const& MeshNodesAlongPolyline::getMesh () const
{
	return _mesh;
}

std::vector<size_t> const& MeshNodesAlongPolyline::getNodeIDs () const
{
	return _msh_node_ids;
}

GeoLib::Polyline const& MeshNodesAlongPolyline::getPolyline () const
{
	return _ply;
}

std::vector<double> const& MeshNodesAlongPolyline::getDistOfProjNodeFromPlyStart() const
{
	return _dist_of_proj_node_from_ply_start;
}
} // end namespace MeshGeoTools
