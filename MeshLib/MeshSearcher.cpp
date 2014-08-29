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

namespace MeshLib
{

std::vector<std::size_t> getConnectedElements(MeshLib::Mesh const& msh, const std::vector<std::size_t> &nodes)
{
	std::vector<std::size_t> connected_elements;
	for (auto* e : msh.getElements()) {
		for (std::size_t j=0; j<e->getNNodes(); j++) {
			if (std::find(nodes.begin(), nodes.end(), e->getNodeIndex(j))!=nodes.end()) {
				connected_elements.push_back(e->getID());
				break;
			}
		}
	}
	return connected_elements;
}

} // end namespace MeshLib

