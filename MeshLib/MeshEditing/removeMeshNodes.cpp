/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of removeMeshNodes.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "removeMeshNodes.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib {

void removeMeshNodes(MeshLib::Mesh &mesh, const std::vector<std::size_t> &del_nodes_idx)
{
	const size_t nDelNodes = del_nodes_idx.size();
	if (nDelNodes == 0)
		return;

	std::vector<MeshLib::Node*>& nodes (mesh._nodes);
	std::vector<MeshLib::Element*>& elements = mesh._elements;
	bool elements_removed (false);

	// delete nodes
	for (std::size_t i = 0; i < nDelNodes; ++i)
	{
		const unsigned idx (del_nodes_idx[i]);
		std::vector<MeshLib::Element*> conn_elems (nodes[idx]->getElements());

		// delete elements connected to these nodes
		for (std::size_t j = 0; j < conn_elems.size(); ++j)
		{
			elements_removed = true;
			auto del_elem (std::find(elements.begin(), elements.end(), conn_elems[j]));
			delete *del_elem;
			*del_elem = nullptr;
		}

		delete nodes[idx];
		nodes[idx] = nullptr;
	}

	// due to element removal neighbourhoods have to be reset and additional nodes
	// might need to be deleted as they are no longer part of any element
	if (elements_removed)
	{
		auto elem_vec_end = std::remove(elements.begin(), elements.end(), nullptr);
		elements.erase(elem_vec_end, elements.end());
		mesh.resetElementsConnectedToNodes();
		for (auto node = nodes.begin(); node != nodes.end(); ++node)
			if ((*node) && (*node)->getNElements() == 0)
			{
				delete *node;
				*node = nullptr;
			}
	}

	auto node_vec_end = std::remove(nodes.begin(), nodes.end(), nullptr);
	nodes.erase(node_vec_end, nodes.end());

	mesh.resetNodeIDs(); // set new node-IDs
}


} // end namespace MeshLib
