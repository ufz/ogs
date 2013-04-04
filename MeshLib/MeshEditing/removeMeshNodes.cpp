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

MeshLib::Mesh* removeMeshNodes(MeshLib::Mesh* mesh, const std::vector<size_t> &nodes)
{
	MeshLib::Mesh* new_mesh (new MeshLib::Mesh(*mesh));

	// delete nodes and their connected elements and replace them with null pointers
	const size_t delNodes = nodes.size();
	std::vector<MeshLib::Node*> mesh_nodes = new_mesh->getNodes();
	for (size_t i = 0; i < delNodes; ++i)
	{
		const MeshLib::Node* node = new_mesh->getNode(i);
		std::vector<MeshLib::Element*> conn_elems = node->getElements();

		for (size_t j = 0; j < conn_elems.size(); ++j)
		{
			delete conn_elems[j];
			conn_elems[j] = NULL;
		}
		delete mesh_nodes[i];
		mesh_nodes[i] = NULL;
	}

	// create map to adjust node indices in element vector
	const size_t nNodes = new_mesh->getNNodes();
	std::vector<int> id_map(nNodes, -1);
	size_t count(0);
	for (size_t i = 0; i < nNodes; ++i)
	{
		if (mesh_nodes[i])
		{
			mesh_nodes[i]->setID(count);
			id_map.push_back(count++);
		}
	}

	// erase null pointers from node- and element vectors
	std::vector<MeshLib::Element*> elements = new_mesh->getElements();
	for (std::vector<MeshLib::Element*>::iterator it = elements.begin(); it != elements.end(); )
	{
		if (*it)
			++it;
		else
			it = elements.erase(it);
	}

	for (std::vector<MeshLib::Node*>::iterator it = mesh_nodes.begin(); it != mesh_nodes.end(); )
	{
		if (*it)
			++it;
		else
			it = mesh_nodes.erase(it);
	}

	return new_mesh;
}

} // end namespace MeshLib
