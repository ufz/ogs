/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MeshCoarsener.cpp
 *
 *  Created on  Aug 3, 2012 by Thomas Fischer
 */

#include "MeshCoarsener.h"

// GeoLib
#include "Grid.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

namespace MeshLib {

MeshCoarsener::MeshCoarsener(Mesh const*const mesh) :
	_orig_mesh (mesh)
{}

MeshCoarsener::~MeshCoarsener()
{}

Mesh* MeshCoarsener::operator()(double min_distance)
{
	// copy mesh nodes, reset mesh node ids and
	// store original mesh node ids in a vector
	std::vector<Node*> const& orig_nodes(_orig_mesh->getNodes());
	const size_t n_nodes(orig_nodes.size());
	std::vector<Node*> nodes(n_nodes);
	std::map<size_t,size_t> orig_ids_map;
	for (size_t k(0); k < n_nodes; k++) {
		nodes[k] = new Node(orig_nodes[k]->getCoords(), k);
		orig_ids_map.insert(std::pair<size_t, size_t>(orig_nodes[k]->getID(), k));
	}

	// init grid
	GeoLib::Grid<Node>* grid(new GeoLib::Grid<Node>(nodes, 64));

	// init id map
	std::vector<size_t> id_map(n_nodes);
	for (size_t k(0); k < n_nodes; k++) {
		id_map[k] = nodes[k]->getID();
	}

	const double sqr_min_distance (min_distance * min_distance);

	// do the work - search nearest nodes
	for (size_t k(0); k < n_nodes; k++) {
		std::vector<std::vector<Node*> const*> node_vecs_intersecting_cube;
		Node const*const node(nodes[k]);
		const size_t node_id(node->getID());
		grid->getVecsOfGridCellsIntersectingCube(node->getCoords(), min_distance, node_vecs_intersecting_cube);

		const size_t n_vecs (node_vecs_intersecting_cube.size());
		for (size_t i(0); i<n_vecs; i++) {
			std::vector<Node*> const* node_vec (node_vecs_intersecting_cube[i]);
			const size_t n_loc_nodes (node_vec->size());
			for (size_t j(0); j<n_loc_nodes; j++) {
				if (node_id < (*node_vec)[j]->getID()) {
					if (MathLib::sqrDist(node->getCoords(), (*node_vec)[j]->getCoords()) < sqr_min_distance) {
						// two nodes are very close to each other
						id_map[j] = k;
					}
				}
			}
		}
	}

	// apply changes to id_map
	for (size_t k(0), cnt(0); k < n_nodes; k++) {
		if (id_map[k] != k) {
			if (id_map[k] != id_map[id_map[k]]) {
				id_map[k] = id_map[id_map[k]];
			}
		} else {
			id_map[k] = cnt++;
		}
	}

	// delete unused nodes
	for (size_t k(0); k < n_nodes; k++) {
		if (id_map[k] != k) {
			delete nodes[k];
			nodes[k] = NULL;
		}
	}

	// remove NULL-ptr from node vector
	for (std::vector<Node*>::iterator it(nodes.begin()); it != nodes.end(); ) {
		if (*it == NULL) {
			it = nodes.erase (it);
		} else {
			it++;
		}
	}

	// copy mesh elements, reset the node pointers
	std::vector<Element*> const& orig_elements(_orig_mesh->getElements());
	const size_t n_elements(orig_elements.size());
	std::vector<Element*> elements(n_elements);
	for (size_t k(0); k < n_elements; k++) {
		elements[k] = orig_elements[k]->clone();
		const size_t n_nodes_element (elements[k]->getNNodes());
		for (size_t i(0); i<n_nodes_element; i++) {
			const size_t orig_node_id (elements[k]->getNode(i)->getID());
			size_t orig_node_pos (std::numeric_limits<size_t>::max());
			std::map<size_t, size_t>::const_iterator it(orig_ids_map.find(orig_node_id));
			if (it == orig_ids_map.end()) {
				std::cerr << "[MeshCoarsener::operator()] could not found mesh node id" << std::endl;
			}
			orig_node_pos = it->second;
			elements[k]->setNode(i, nodes[id_map[orig_node_pos]]);
		}
	}

	return new Mesh ("test", nodes, elements);
}

} // end namespace MeshLib
