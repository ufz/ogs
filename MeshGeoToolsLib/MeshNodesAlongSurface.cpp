/**
 * \author Norihiro Watanabe
 * \date   2014-03-14
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshNodesAlongSurface.h"

#include <algorithm>

#include "quicksort.h"
#include "MathTools.h"


namespace MeshGeoToolsLib
{

MeshNodesAlongSurface::MeshNodesAlongSurface(
		MeshLib::Mesh const& mesh,
		GeoLib::Surface const& sfc) :
	_mesh(mesh), _sfc(sfc)
{
	auto& mesh_nodes = _mesh.getNodes();
	const std::size_t n_nodes (mesh_nodes.size());
	// loop over all nodes
	for (size_t i = 0; i < n_nodes; i++) {
		auto* node = mesh_nodes[i];
		if (!sfc.isPntInBoundingVolume(node->getCoords()))
			continue;
		if (sfc.isPntInSfc(node->getCoords())) {
			_msh_node_ids.push_back(node->getID());
		}
	}
}

MeshLib::Mesh const& MeshNodesAlongSurface::getMesh () const
{
	return _mesh;
}

std::vector<std::size_t> const& MeshNodesAlongSurface::getNodeIDs () const
{
	return _msh_node_ids;
}

GeoLib::Surface const& MeshNodesAlongSurface::getSurface () const
{
	return _sfc;
}

} // end namespace MeshGeoTools
