/**
 * \file   AddLayerToMesh.cpp
 * \date   2016-01-18
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
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEditing/FlipElements.h"

namespace MeshLib
{

MeshLib::Element* extrudeElement(std::vector<MeshLib::Node*> const& subsfc_nodes,
	MeshLib::Element const& sfc_elem,
	std::map<std::size_t, std::size_t> const& subsfc_sfc_id_map)
{
	if (sfc_elem.getDimension() > 2)
		return nullptr;

	const unsigned nElemNodes(sfc_elem.getNBaseNodes());
	MeshLib::Node** new_nodes = new MeshLib::Node*[2*nElemNodes];

	for (unsigned j=0; j<nElemNodes; ++j)
	{
		new_nodes[j] = subsfc_nodes[sfc_elem.getNode(j)->getID()];
		std::size_t new_idx = (nElemNodes==2) ? (3-j) : (nElemNodes+j);
		new_nodes[new_idx] =
		    subsfc_nodes[subsfc_sfc_id_map.at(sfc_elem.getNode(j)->getID())];
	}
	
	if (sfc_elem.getGeomType() == MeshLib::MeshElemType::LINE)
		return new MeshLib::Quad(new_nodes);
	if (sfc_elem.getGeomType() == MeshLib::MeshElemType::TRIANGLE)
		return new MeshLib::Prism(new_nodes);
	if (sfc_elem.getGeomType() == MeshLib::MeshElemType::QUAD)
		return new MeshLib::Hex(new_nodes);

	return nullptr;
}

MeshLib::Mesh* addTopLayerToMesh(MeshLib::Mesh const& mesh,
	double thickness,
	std::string const& name)
{
	return addLayerToMesh(mesh, thickness, name, true);
}

MeshLib::Mesh* addBottomLayerToMesh(MeshLib::Mesh const& mesh,
	double thickness,
	std::string const& name)
{
	return addLayerToMesh(mesh, thickness, name, false);
}

MeshLib::Mesh* addLayerToMesh(MeshLib::Mesh const& mesh, double thickness,
	std::string const& name,
	bool on_top)
{
	INFO("Extracting top surface of mesh \"%s\" ... ", mesh.getName().c_str());
	int const flag = (on_top) ? -1 : 1;
	const MathLib::Vector3 dir(0, 0, flag);
	double const angle(90);
	std::unique_ptr<MeshLib::Mesh> sfc_mesh (nullptr);
	
	if (mesh.getDimension() == 3)
		sfc_mesh.reset(MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, dir, angle, true));
	else
		sfc_mesh = (on_top) ? std::unique_ptr<MeshLib::Mesh>(new MeshLib::Mesh(mesh)) :
		                      std::unique_ptr<MeshLib::Mesh>(MeshLib::createFlippedMesh(mesh));
	INFO("done.");

	// *** add new surface nodes
	std::vector<MeshLib::Node*> subsfc_nodes =
	    MeshLib::copyNodeVector(mesh.getNodes());
	std::vector<MeshLib::Element*> subsfc_elements =
		MeshLib::copyElementVector(mesh.getElements(), subsfc_nodes);

	std::size_t const n_subsfc_nodes(subsfc_nodes.size());

	std::vector<MeshLib::Node*> const& sfc_nodes(sfc_mesh->getNodes());
	std::size_t const n_sfc_nodes(sfc_nodes.size());

	// *** copy sfc nodes to subsfc mesh node
	std::map<std::size_t, std::size_t> subsfc_sfc_id_map;
	for (std::size_t k(0); k<n_sfc_nodes; ++k) {
		std::size_t const subsfc_id(sfc_nodes[k]->getID());
		std::size_t const sfc_id(k+n_subsfc_nodes);
		subsfc_sfc_id_map.insert(std::make_pair(subsfc_id, sfc_id));
		MeshLib::Node const& node (*sfc_nodes[k]);
		subsfc_nodes.push_back(new MeshLib::Node(
		    node[0], node[1], node[2] - (flag * thickness), sfc_id));
	}

	// *** insert new top layer elements into subsfc_mesh
	std::vector<MeshLib::Element*> const& sfc_elements(sfc_mesh->getElements());
	std::size_t const n_sfc_elements(sfc_elements.size());
	for (std::size_t k(0); k<n_sfc_elements; ++k)
		subsfc_elements.push_back(
			extrudeElement(subsfc_nodes, *sfc_elements[k], subsfc_sfc_id_map)
		);

	auto new_mesh = new MeshLib::Mesh(name, subsfc_nodes, subsfc_elements);

	boost::optional<MeshLib::PropertyVector<int> const&> opt_materials(
		mesh.getProperties().getPropertyVector<int>("MaterialIDs")
	);

	if (opt_materials) {
		boost::optional<PropertyVector<int> &> new_opt_materials(
		new_mesh->getProperties().createNewPropertyVector<int>("MaterialIDs",
			MeshLib::MeshItemType::Cell, 1));
		if (!new_opt_materials) {
			ERR("Can not set material properties for new layer");
		} else {
			unsigned new_layer_id(*(std::max_element(
				opt_materials->cbegin(), opt_materials->cend()))+1);
			std::copy(opt_materials->cbegin(), opt_materials->cend(),
			          new_opt_materials->begin());
			auto const n_new_props(subsfc_elements.size()-mesh.getNElements());
			std::fill_n(new_opt_materials->end(), n_new_props, new_layer_id);
		}
	} else {
		ERR(
		    "Could not copy the property \"MaterialIDs\" since the original "
		    "mesh does not contain such a property.");
	}

	return new_mesh;
}

} // namespace MeshLib
