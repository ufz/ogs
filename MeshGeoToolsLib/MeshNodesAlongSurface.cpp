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
		std::vector<MeshLib::Node*> const& mesh_nodes,
		GeoLib::Surface const& sfc) :
	_sfc(sfc)
{
	for (auto itr=mesh_nodes.cbegin(); itr!=mesh_nodes.cend(); ++itr) {
		auto node = *itr;
		if (!sfc.isPntInBoundingVolume(node->getCoords()))
			continue;
		if (sfc.isPntInSfc(node->getCoords())) {
			_msh_node_ids.push_back(node->getID());
		}
	}
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
