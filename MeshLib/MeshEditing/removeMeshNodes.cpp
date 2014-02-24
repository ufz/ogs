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
	// delete nodes and their connected elements and replace them with null pointers
	const size_t nDelNodes = del_nodes_idx.size();
	std::vector<MeshLib::Node*>& mesh_nodes = mesh.getNodesVector();
	
	for (std::size_t i = 0; i < nDelNodes; ++i)
	{
		const unsigned idx (del_nodes_idx[i]);
		std::vector<MeshLib::Element*> conn_elems (mesh_nodes[idx]->getElements());

		for (std::size_t j = 0; j < conn_elems.size(); ++j)
		{
			delete conn_elems[j];
			conn_elems[j] = nullptr;
		}
		delete mesh_nodes[idx];
		mesh_nodes[idx] = nullptr;
	}

	// erase null pointers from node- and element vectors
	std::vector<MeshLib::Element*> elements = mesh.getElements();
	auto elem_vec_end = std::remove(elements.begin(), elements.end(), nullptr);
	elements.erase(elem_vec_end, elements.end());

	auto node_vec_end = std::remove(mesh_nodes.begin(), mesh_nodes.end(), nullptr);
	mesh_nodes.erase(node_vec_end, mesh_nodes.end());
	
	if (nDelNodes)
		mesh.resetNodeIDs(); // after removing nodes set new node-IDs
}


} // end namespace MeshLib
