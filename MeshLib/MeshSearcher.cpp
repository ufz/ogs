/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

std::vector<std::size_t> getConnectedNodeIDs(const std::vector<MeshLib::Element*> &elements)
{
	std::vector<std::size_t> connected_nodes;
	std::for_each(elements.begin(), elements.end(),
		[&](MeshLib::Element* e)
		{
			for (unsigned i=0; i<e->getNNodes(); i++)
				connected_nodes.push_back(e->getNodeIndex(i));
		});
	std::sort(connected_nodes.begin(), connected_nodes.end());
	auto it = std::unique(connected_nodes.begin(), connected_nodes.end());
	connected_nodes.resize(std::distance(connected_nodes.begin(),it));
	return connected_nodes;
}

} // end namespace MeshLib

