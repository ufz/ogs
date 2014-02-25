/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of removeMeshNodes.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "removeMeshNodes.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

namespace MeshLib {

void removeMeshNodes(MeshLib::Mesh &mesh, const std::vector<std::size_t> &del_nodes_idx)
{
	const size_t nDelNodes = del_nodes_idx.size();
	if (nDelNodes == 0)
		return;

	// delete nodes and their connected elements and replace them with null pointers
	std::vector<MeshLib::Node*>& nodes (mesh._nodes);
	std::vector<MeshLib::Element*>& elements = mesh._elements;

	for (std::size_t i = 0; i < nDelNodes; ++i)
	{
		const unsigned idx (del_nodes_idx[i]);
		std::vector<MeshLib::Element*> conn_elems (nodes[idx]->getElements());
		
		for (std::size_t j = 0; j < conn_elems.size(); ++j)
		{
			auto del_elem (std::find(elements.begin(), elements.end(), conn_elems[j]));
			delete *del_elem;
			*del_elem = nullptr;
		}
		
		delete nodes[idx];
		nodes[idx] = nullptr;
	}
	
	auto elem_vec_end = std::remove(elements.begin(), elements.end(), nullptr);
	elements.erase(elem_vec_end, elements.end());
	auto node_vec_end = std::remove(nodes.begin(), nodes.end(), nullptr);
	nodes.erase(node_vec_end, nodes.end());

	mesh.resetNodeIDs(); // after removing nodes set new node-IDs
}


} // end namespace MeshLib
