/**
 * \file
 * \author Thomas Fischer
 * \date   Oct 12, 2012
 * \brief  Implementation of the Mesh2MeshPropertyInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>
#include <fstream>
#include <boost/optional.hpp>

#include "Mesh2MeshPropertyInterpolation.h"

#include "logog/include/logog.hpp"

#include "GeoLib/AABB.h"
#include "GeoLib/Grid.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

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
		ERR("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() source mesh to small.");
		ERR("src_aabb: %f, %f, %f | %f, %f, %f", src_aabb.getMinPoint()[0], src_aabb.getMinPoint()[1], src_aabb.getMinPoint()[2], src_aabb.getMaxPoint()[0], src_aabb.getMaxPoint()[1], src_aabb.getMaxPoint()[2]);
		ERR("dest_aabb: %f, %f, %f | %f, %f, %f", dest_aabb.getMinPoint()[0], dest_aabb.getMinPoint()[1], dest_aabb.getMinPoint()[2], dest_aabb.getMaxPoint()[0], dest_aabb.getMaxPoint()[1], dest_aabb.getMaxPoint()[2]);
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

	GeoLib::Grid<MeshLib::Node> src_grid(src_nodes.begin(), src_nodes.end(), 64);

	auto materialIds = dest_mesh->getProperties().getPropertyVector<int>("MaterialIDs");
	if (!materialIds)
		materialIds = dest_mesh->getProperties().createNewPropertyVector<int>(
			"MaterialIDs", MeshLib::MeshItemType::Cell, 1);
	if (!materialIds)
	{
		ERR("Could not create PropertyVector for MaterialIDs in Mesh.");
		return;
	}

	std::vector<MeshLib::Element*> const& dest_elements(dest_mesh->getElements());
	const std::size_t n_dest_elements(dest_elements.size());
	for (std::size_t k(0); k<n_dest_elements; k++)
	{
		// compute axis aligned bounding box around the current element
		const GeoLib::AABB<MeshLib::Node> elem_aabb(dest_elements[k]->getNodes(), dest_elements[k]->getNodes()+dest_elements[k]->getNBaseNodes());

		// request "interesting" nodes from grid
		std::vector<std::vector<MeshLib::Node*> const*> nodes;
		src_grid.getPntVecsOfGridCellsIntersectingCuboid(elem_aabb.getMinPoint(), elem_aabb.getMaxPoint(), nodes);

		std::size_t cnt(0);
		dest_properties[k] = 0.0;

		for (std::size_t i(0); i<nodes.size(); ++i) {
			std::vector<MeshLib::Node*> const* i_th_vec(nodes[i]);
			const std::size_t n_nodes_in_vec(i_th_vec->size());
			for (std::size_t j(0); j<n_nodes_in_vec; j++) {
				MeshLib::Node const*const j_th_node((*i_th_vec)[j]);
				if (elem_aabb.containsPoint(*j_th_node)) {
					if (dest_elements[k]->isPntInElement(*j_th_node, 30)) {
						dest_properties[k] += interpolated_src_node_properties[(*i_th_vec)[j]->getID()];
						cnt++;
					}
				}
			}
		}

		dest_properties[k] /= cnt;
		materialIds->push_back(k);

		if (cnt == 0) {
			std::string base_name("DebugData/Element-");
			base_name += std::to_string(k);

			std::string aabb_fname(base_name + "-aabb.gli");
			std::ofstream out_aabb(aabb_fname.c_str());
			out_aabb << "#POINTS" << "\n";
			out_aabb << "0 " << elem_aabb.getMinPoint() << "\n";
			out_aabb << "1 " << elem_aabb.getMaxPoint() << "\n";
			out_aabb << "#STOP" << "\n";
			out_aabb.close();


			std::string source_fname(base_name + "-SourceNodes.gli");
			std::ofstream out_src(source_fname.c_str());
			out_src << "#POINTS" << "\n";
			std::size_t nodes_cnt(0);
			for (std::size_t i(0); i<nodes.size(); ++i) {
				std::vector<MeshLib::Node*> const* i_th_vec(nodes[i]);
				const std::size_t n_nodes_in_vec(i_th_vec->size());
				for (std::size_t j(0); j<n_nodes_in_vec; j++) {
					MeshLib::Node const*const j_th_node((*i_th_vec)[j]);
					out_src << nodes_cnt << " " << *(dynamic_cast<GeoLib::Point const*>(j_th_node)) << "\n";
					nodes_cnt++;
				}
			}
			out_src << "#STOP" << "\n";
			out_src.close();
			std::cerr << "no source nodes in dest element " << k << "\n";
		}
	}
}

void Mesh2MeshPropertyInterpolation::interpolateElementPropertiesToNodeProperties(std::vector<double> &interpolated_node_properties) const
{
	auto materialIds = _src_mesh->getProperties().getPropertyVector<int>("MaterialIDs");
	if (!materialIds)
		return;

	std::vector<MeshLib::Node*> const& src_nodes(_src_mesh->getNodes());
	const std::size_t n_src_nodes(src_nodes.size());
	for (std::size_t k(0); k<n_src_nodes; k++) {
		const std::size_t n_con_elems (src_nodes[k]->getNElements());
		interpolated_node_properties[k] = (*_src_properties)[(*materialIds)[src_nodes[k]->getElement(0)->getID()]];
		for (std::size_t j(1); j<n_con_elems; j++) {
			interpolated_node_properties[k] += (*_src_properties)[(*materialIds)[src_nodes[k]->getElement(j)->getID()]];
		}
		interpolated_node_properties[k] /= n_con_elems;
	}
}

} // end namespace MeshLib
