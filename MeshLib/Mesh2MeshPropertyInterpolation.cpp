/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file Mesh2MeshPropertyInterpolation.cpp
 *
 *  Created on  Oct 12, 2012 by Thomas Fischer
 */

#include <vector>
#include <fstream>

#include "Mesh2MeshPropertyInterpolation.h"

// BaseLib
#include "logog.hpp"

// GeoLib
#include "AABB.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Face.h"

namespace MeshLib {

Mesh2MeshPropertyInterpolation::Mesh2MeshPropertyInterpolation(Mesh const*const src_mesh, std::vector<double> const*const src_properties) :
	_src_mesh(src_mesh), _src_properties(src_properties)
{}

Mesh2MeshPropertyInterpolation::~Mesh2MeshPropertyInterpolation()
{}

bool Mesh2MeshPropertyInterpolation::setPropertiesForMesh(Mesh *dest_mesh, std::vector<double>& dest_properties) const
{
	if (_src_mesh->getDimension() != dest_mesh->getDimension()) {
		ERR ("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() dimension of source (dim = %d) and destination (dim = %d) mesh does not match.", _src_mesh->getDimension(), dest_mesh->getDimension());
		return false;
	}

	if (_src_mesh->getDimension() != 2) {
		WARN ("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() implemented only for 2D case at the moment.");
		return false;
	}

	GeoLib::AABB<MeshLib::Node> src_aabb(_src_mesh->getNodes().begin(), _src_mesh->getNodes().end());
	GeoLib::AABB<MeshLib::Node> dest_aabb(dest_mesh->getNodes().begin(), dest_mesh->getNodes().end());
	if (!src_aabb.containsAABB(dest_aabb)) {
		ERR ("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() source mesh to small.");
		return false;
	}

	interpolatePropertiesForMesh(dest_mesh, dest_properties);

	return true;
}

void Mesh2MeshPropertyInterpolation::interpolatePropertiesForMesh(Mesh *dest_mesh, std::vector<double>& dest_properties) const
{
	// carry over property information from source elements to source nodes
	std::vector<double> interpolated_src_node_properties(_src_mesh->getNNodes());
	interpolateElementPropertiesToNodeProperties(interpolated_src_node_properties);

	// looping over the destination elements and calculate properties
	// from interpolated_src_node_properties
	std::vector<MeshLib::Node*> const& src_nodes(_src_mesh->getNodes());
	const size_t n_src_nodes(src_nodes.size());
	std::vector<MeshLib::Element*> const& dest_elements(dest_mesh->getElements());
	const size_t n_dest_elements(dest_elements.size());
	for (size_t k(0); k<n_dest_elements; k++) {
		size_t cnt(0);
		dest_properties[k] = 0.0;
		for (size_t j(0); j<n_src_nodes; j++) {
			if (dynamic_cast<MeshLib::Face*>(dest_elements[k])->isPntInside(* dynamic_cast<GeoLib::Point const*>(src_nodes[j]))) {
				dest_properties[k] += interpolated_src_node_properties[j];
				cnt++;
			}
		}

		dest_properties[k] /= cnt;
		dest_elements[k]->setValue(k);
	}
}

void Mesh2MeshPropertyInterpolation::interpolateElementPropertiesToNodeProperties(std::vector<double> &interpolated_node_properties) const
{
	std::vector<MeshLib::Node*> const& src_nodes(_src_mesh->getNodes());
	const size_t n_src_nodes(src_nodes.size());

	for (size_t k(0); k<n_src_nodes; k++) {
		const size_t n_con_elems (src_nodes[k]->getNElements());
		interpolated_node_properties[k] = (*_src_properties)[(src_nodes[k]->getElement(0))->getValue()];
		for (size_t j(1); j<n_con_elems; j++) {
			interpolated_node_properties[k] += (*_src_properties)[(src_nodes[k]->getElement(j))->getValue()];
		}
		interpolated_node_properties[k] /= n_con_elems;
	}
}

} // end namespace MeshLib
