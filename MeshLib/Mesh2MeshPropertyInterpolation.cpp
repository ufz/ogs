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
#include "AxisAlignedBoundingBox.h"

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

	GeoLib::AABB src_aabb(_src_mesh->getNodes().begin(), _src_mesh->getNodes().end());
	GeoLib::AABB dest_aabb(dest_mesh->getNodes().begin(), dest_mesh->getNodes().end());
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

	// carry over properties from source elements to destination nodes
	// stepping over source nodes of the appropriate source element
	std::vector<MeshLib::Node*> const& dest_nodes(dest_mesh->getNodes());
	const size_t n_dest_nodes(dest_nodes.size());
	std::vector<double> dest_node_properties(n_dest_nodes);
	std::vector<MeshLib::Element*> const& src_elements(_src_mesh->getElements());
	const size_t n_src_elements(src_elements.size());
	for (size_t k(0); k<n_dest_nodes; k++) {
		// search the source element the destination node is in
		for (size_t j(0); j<n_src_elements; j++) {
			if (dynamic_cast<MeshLib::Face*>(src_elements[j])->isPntInside(* dynamic_cast<GeoLib::Point const*>(dest_nodes[k]))) {
				const size_t n_nodes_src_element(src_elements[j]->getNNodes());
				dest_node_properties[k] = interpolated_src_node_properties[(src_elements[j]->getNode(0))->getID()];
				for (size_t i(1); i<n_nodes_src_element; i++) {
					dest_node_properties[k] += interpolated_src_node_properties[(src_elements[j]->getNode(i))->getID()];
				}
			}
		}
	}

	// looping over the destination elements and interpolate properties
	// from dest_node_properties
	std::vector<MeshLib::Element*> const& dest_elements(dest_mesh->getElements());
	const size_t n_dest_elements(dest_elements.size());
	for (size_t k(0); k<n_dest_elements; k++) {
		const size_t n_nodes_dest_element(dest_elements[k]->getNNodes());
		dest_properties[k] = dest_node_properties[dest_elements[k]->getNode(0)->getID()];
		for (size_t j(1); k<n_nodes_dest_element; k++) {
			dest_properties[k] += dest_node_properties[dest_elements[k]->getNode(j)->getID()];
		}
		dest_properties[k] /= n_nodes_dest_element;
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

//void Mesh2MeshPropertyInterpolation::interpolatePropertiesForMeshUsingArea(Mesh *dest_mesh, std::vector<double>& dest_properties) const
//{
//	// interpolate properties for all all elements of dest_mesh
//	std::vector<MeshLib::Element*> const* dest_elements(dest_mesh->getElements());
//	const size_t n_dest_elements(dest_elements.size());
//	for (size_t k(0); k<n_dest_elements; k++) {
//		std::vector<double> interpolation_weights;
//		std::vector<double> interpolation_values;
//		// find all elements of _src_mesh that have a non empty intersection and compute the (percentaged) part of volume
//		{
//			const double content_dest_elem(dest_elements[k]->getContent()); // needed for stop criterion
//			double sum_of_interpolation_weights(0.0);
//			std::vector<MeshLib::Element*> const* src_elements(_src_mesh->getElements());
//			const size_t n_src_elements (src_elements.size());
//			for (size_t j(0); j<n_src_elements && sum_of_interpolation_weights <= 0.999; j++) {
//				double act_interpolation_weight(getIntersectingContent(dest_elements[k], src_elements[j]) / content_dest_elem);
//				if (act_interpolation_weight > 0) {
//					sum_of_interpolation_weights += act_interpolation_weight;
//					interpolation_weights.push_back(act_interpolation_weight);
//					interpolation_values.push_back(_src_properties[src_elements[j]->getValue()]);
//				}
//			}
//		}
//
//		// create the new property out of the interpolation of properties of the above results
//		const size_t n_interpolation_info(interpolation_weights.size());
//		double new_prop (0.0);
//		for (size_t i(0); i<n_interpolation_info;i++) {
//			new_prop += interpolation_weights[i] * interpolation_values[i];
//		}
//		dest_properties.push_back(new_prop);
//
//		// set index of element to the index the new property resides in the vector
//		const_cast<MeshLib::Element*>(dest_elements[k])->setValue(dest_properties.size()-1);
//	}
//}

//double Mesh2MeshPropertyInterpolation::getIntersectingContent(Element const*const elem0, Element const*const elem1) const
//{
//	const unsigned n_edges_0(elem0->getNEdges());
//	const unsigned n_edges_1(elem1->getNEdges());
//
//	for (unsigned i(0); i<n_edges_0; i++) {
//		for (unsigned j(0); k<n_edges_1; j++) {
//
//		}
//	}
//}

} // end namespace MeshLib
