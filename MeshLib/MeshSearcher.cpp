/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshSearcher.h"

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

#include <algorithm>

namespace MeshLib
{

std::vector<std::size_t> getConnectedElementIDs(MeshLib::Mesh const& msh, const std::vector<std::size_t> &nodes)
{
	std::vector<std::size_t> connected_elements;
	std::for_each(nodes.begin(), nodes.end(),
		[&](std::size_t node_id)
		{
			for (auto* e : msh.getNode(node_id)->getElements()) {
				connected_elements.push_back(e->getID());
			}
		});
	std::sort(connected_elements.begin(), connected_elements.end());
	auto it = std::unique(connected_elements.begin(), connected_elements.end());
	connected_elements.resize(std::distance(connected_elements.begin(),it));
	return connected_elements;
}

} // end namespace MeshLib

