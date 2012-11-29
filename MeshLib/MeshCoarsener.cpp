/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file MeshCoarsener.cpp
 *
 *  Created on  Aug 3, 2012 by Thomas Fischer
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

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
	GeoLib::Grid<Node>* grid(new GeoLib::Grid<Node>(nodes.begin(), nodes.end(), 64));

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
		const size_t node_id_k(node->getID());
		grid->getVecsOfGridCellsIntersectingCube(node->getCoords(), min_distance, node_vecs_intersecting_cube);

		const size_t n_vecs (node_vecs_intersecting_cube.size());
		for (size_t i(0); i<n_vecs; i++) {
			std::vector<Node*> const* node_vec (node_vecs_intersecting_cube[i]);
			const size_t n_loc_nodes (node_vec->size());
			for (size_t j(0); j<n_loc_nodes; j++) {
				Node const*const test_node((*node_vec)[j]);
				const size_t test_node_id (test_node->getID());
				if (node_id_k < test_node_id) {
					if (MathLib::sqrDist(node->getCoords(), test_node->getCoords()) < sqr_min_distance) {
						// two nodes are very close to each other
						id_map[test_node_id] = node_id_k;
#ifndef NDEBUG
						INFO ("distance of nodes with ids %d and %d is %f", node_id_k, test_node_id, sqrt(MathLib::sqrDist(node->getCoords(), test_node->getCoords())));
#endif
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
	for (size_t k(0), cnt(0); k < n_nodes; k++) {
		if (id_map[k] != cnt) {
			delete nodes[k];
			nodes[k] = NULL;
		} else {
			cnt++;
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

	// reset mesh node ids
	const size_t new_n_nodes(nodes.size());
	for (size_t k(0); k < new_n_nodes; k++) {
		nodes[k]->setID(k);
	}

	// copy mesh elements, reset the node pointers
	std::vector<Element*> const& orig_elements(_orig_mesh->getElements());
	const size_t n_elements(orig_elements.size());
	std::vector<Element*> elements(n_elements);
	std::vector<size_t> mapped_node_ids_of_element;
	for (size_t k(0), cnt(0); k < n_elements; k++) {
		Element const*const kth_orig_elem(orig_elements[k]);
		const size_t n_nodes_element (kth_orig_elem->getNNodes());
		mapped_node_ids_of_element.clear();
		for (size_t i(0); i<n_nodes_element; i++) {
			const size_t orig_node_id (kth_orig_elem->getNode(i)->getID());
			std::map<size_t, size_t>::const_iterator it(orig_ids_map.find(orig_node_id));
			if (it == orig_ids_map.end()) {
				std::cerr << "[MeshCoarsener::operator()] could not found mesh node id" << std::endl;
			} else {
				mapped_node_ids_of_element.push_back(id_map[it->second]);
			}
		}

		// check if nodes of the element are collapsed
		bool not_collapsed (true);
		for (size_t i(0); i<n_nodes_element-1 && not_collapsed; i++) {
			const size_t id_i(mapped_node_ids_of_element[i]);
			for (size_t j(i+1); j<n_nodes_element && not_collapsed; j++) {
				if (id_i == mapped_node_ids_of_element[j]) {
					not_collapsed = false;
				}
			}
		}
		if (! not_collapsed) {
			Element* elem (kth_orig_elem->clone());
			if (elem != NULL) {
				for (size_t i(0); i<n_nodes_element; i++) {
					elem->setNode(i, nodes[mapped_node_ids_of_element[i]]);
				}
				Element* revised_elem(elem->reviseElement());
				elements[cnt] = revised_elem;
				delete elem;
				cnt++;
			}
		} else {
			elements[cnt] = kth_orig_elem->clone();
			for (size_t i(0); i<n_nodes_element; i++) {
				elements[cnt]->setNode(i, nodes[mapped_node_ids_of_element[i]]);
			}
			cnt++;
		}
	}

	for (std::vector<Element*>::iterator it(elements.begin()); it != elements.end(); ) {
		if (*it == NULL) {
			it = elements.erase(it);
		} else {
			it++;
		}
	}

	return new Mesh (_orig_mesh->getName() + "Collapsed", nodes, elements);
}

} // end namespace MeshLib
