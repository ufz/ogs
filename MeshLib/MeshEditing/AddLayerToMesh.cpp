/**
 * \file   AddLayerToMesh.cpp
 * \author Karsten Rink
 * \date   2014-03-25
 * \brief  Implementation of AddLayerToMesh class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AddLayerToMesh.h"

#include <vector>
#include <map>
#include <memory>

#include "logog/include/logog.hpp"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/MeshSurfaceExtraction.h"

namespace MeshLib
{

MeshLib::Prism* extrudeElement(std::vector<MeshLib::Node*> const& subsfc_nodes,
	MeshLib::Tri const*const sfc_elem,
	std::map<std::size_t, std::size_t> const& subsfc_sfc_id_map)
{
	std::array<MeshLib::Node*, 6> prism_nodes;
	prism_nodes[0] = subsfc_nodes[sfc_elem->getNode(0)->getID()];
	prism_nodes[1] = subsfc_nodes[sfc_elem->getNode(1)->getID()];
	prism_nodes[2] = subsfc_nodes[sfc_elem->getNode(2)->getID()];
	prism_nodes[3] = subsfc_nodes[
			subsfc_sfc_id_map.at(sfc_elem->getNode(0)->getID())];
	prism_nodes[4] = subsfc_nodes[
			subsfc_sfc_id_map.at(sfc_elem->getNode(1)->getID())];
	prism_nodes[5] = subsfc_nodes[
				subsfc_sfc_id_map.at(sfc_elem->getNode(2)->getID())];
	return new MeshLib::Prism(prism_nodes);
}

MeshLib::Hex* extrudeElement(std::vector<MeshLib::Node*> const& subsfc_nodes,
	MeshLib::Quad const*const sfc_elem,
	std::map<std::size_t, std::size_t> const& subsfc_sfc_id_map)
{
	std::array<MeshLib::Node*, 8> hex_nodes;
	hex_nodes[0] = subsfc_nodes[sfc_elem->getNode(0)->getID()];
	hex_nodes[1] = subsfc_nodes[sfc_elem->getNode(1)->getID()];
	hex_nodes[2] = subsfc_nodes[sfc_elem->getNode(2)->getID()];
	hex_nodes[3] = subsfc_nodes[sfc_elem->getNode(3)->getID()];
	hex_nodes[4] = subsfc_nodes[
		subsfc_sfc_id_map.at(sfc_elem->getNode(0)->getID())];
	hex_nodes[5] = subsfc_nodes[
		subsfc_sfc_id_map.at(sfc_elem->getNode(1)->getID())];
	hex_nodes[6] = subsfc_nodes[
		subsfc_sfc_id_map.at(sfc_elem->getNode(2)->getID())];
	hex_nodes[7] = subsfc_nodes[
		subsfc_sfc_id_map.at(sfc_elem->getNode(3)->getID())];
	return new MeshLib::Hex(hex_nodes);
}

MeshLib::Quad* extrudeElement(std::vector<MeshLib::Node*> const& subsfc_nodes,
	MeshLib::Line const*const sfc_elem,
	std::map<std::size_t, std::size_t> const& subsfc_sfc_id_map)
{
	std::array<MeshLib::Node*, 4> quad_nodes;
	quad_nodes[0] = subsfc_nodes[sfc_elem->getNode(1)->getID()];
	quad_nodes[1] = subsfc_nodes[sfc_elem->getNode(0)->getID()];
	quad_nodes[2] = subsfc_nodes[
		subsfc_sfc_id_map.at(sfc_elem->getNode(0)->getID())];
	quad_nodes[3] = subsfc_nodes[
		subsfc_sfc_id_map.at(sfc_elem->getNode(1)->getID())];
	return new MeshLib::Quad(quad_nodes);
}

MeshLib::Mesh* addTopLayerToMesh(MeshLib::Mesh const& mesh, double thickness)
{
	INFO("Extracting top surface of mesh \"%s\" ... ", mesh.getName().c_str());
	const MathLib::Vector3 dir(0,0,-1);
	double const angle(90);
	std::unique_ptr<MeshLib::Mesh> sfc_mesh(
		MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, dir, angle, true)
	);
	INFO("done.");

	MeshLib::Mesh* subsfc_mesh = new MeshLib::Mesh (mesh);

	// *** add new top surface nodes
	std::vector<MeshLib::Node*> & subsfc_nodes(
		const_cast<std::vector<MeshLib::Node*> &>(subsfc_mesh->getNodes())
	);
	std::size_t const n_subsfc_nodes(subsfc_nodes.size());

	std::vector<MeshLib::Node*> const& sfc_nodes(sfc_mesh->getNodes());
	std::size_t const n_sfc_nodes(sfc_nodes.size());

	// *** copy sfc nodes to subsfc mesh node
	std::map<std::size_t, std::size_t> subsfc_sfc_id_map;
	for (std::size_t k(0); k<n_sfc_nodes; ++k) {
		subsfc_nodes.push_back(new MeshLib::Node(*sfc_nodes[k]));
		std::size_t const subsfc_id(sfc_nodes[k]->getID());
		std::size_t const sfc_id(k+n_subsfc_nodes);
		subsfc_sfc_id_map.insert(std::make_pair(subsfc_id, sfc_id));
		(*(subsfc_nodes.back()))[2] += thickness;
	}
	subsfc_mesh->resetNodeIDs();

	// *** insert new top layer elements into subsfc_mesh
	std::vector<MeshLib::Element*> & subsfc_elements(
		const_cast<std::vector<MeshLib::Element*> &>(subsfc_mesh->getElements())
	);
	std::size_t orig_size(subsfc_elements.size());
	std::vector<MeshLib::Element*> const& sfc_elements(sfc_mesh->getElements());
	std::size_t const n_sfc_elements(sfc_elements.size());
	for (std::size_t k(0); k<n_sfc_elements; ++k) {
		MeshLib::Element const*const sfc_elem(sfc_elements[k]);
		if (sfc_elem->getGeomType() == MeshLib::MeshElemType::TRIANGLE) {
			// add a prism
			subsfc_elements.push_back(extrudeElement(subsfc_nodes,
				static_cast<MeshLib::Tri const*const>(sfc_elem),
				subsfc_sfc_id_map));
		} else {
			if (sfc_elements[k]->getGeomType() == MeshLib::MeshElemType::QUAD) {
				// add a hexahedron
				subsfc_elements.push_back(extrudeElement(subsfc_nodes,
					static_cast<MeshLib::Quad const*const>(sfc_elem),
					subsfc_sfc_id_map));
			} else {
				if (sfc_elements[k]->getGeomType() == MeshLib::MeshElemType::LINE) {
					// add a quad
					subsfc_elements.push_back(extrudeElement(subsfc_nodes,
						static_cast<MeshLib::Line const*const>(sfc_elem),
						subsfc_sfc_id_map));
				}
			}
		}
	}

	boost::optional<MeshLib::PropertyVector<int> &> opt_materials(
		subsfc_mesh->getProperties().getPropertyVector<int>("MaterialIDs")
	);
	if (!opt_materials) {
		ERR("Can not set material properties for new layer");
	} else {
		MeshLib::PropertyVector<int> & materials(opt_materials.get());
		unsigned layer_id(*(std::max_element(
			materials.cbegin(), materials.cend()))+1);
		while (orig_size<subsfc_elements.size()) {
			materials.push_back(layer_id);
			orig_size++;
		}
	}

	return subsfc_mesh;
}

} // namespace MeshLib
