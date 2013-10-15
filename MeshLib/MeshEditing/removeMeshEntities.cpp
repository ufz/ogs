/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of removeMeshEntities.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "removeMeshEntities.h"
#include "Mesh.h"
#include "Elements/Element.h"
#include "AABB.h"

#include "logog/include/logog.hpp"

namespace MeshLib {

MeshElementRemoval::MeshElementRemoval(const MeshLib::Mesh &mesh)
	: _mesh(mesh)
{
}

MeshElementRemoval::~MeshElementRemoval()
{
}

MeshLib::Mesh* MeshElementRemoval::removeMeshElements() const
{
	INFO("Removing total %d elements...", _marked_elements.size());
	std::vector<MeshLib::Element*> tmp_elems = excludeElements(_mesh.getElements(), _marked_elements);
	INFO("%d elements remain in mesh.", tmp_elems.size());
	std::vector<MeshLib::Node*> new_nodes;
	std::vector<MeshLib::Element*> new_elems;
	copyNodesElements(_mesh.getNodes(), tmp_elems, new_nodes, new_elems);

	// create a new mesh object. Unsued nodes are removed during construction
	return new MeshLib::Mesh(_mesh.getName(), new_nodes, new_elems);
}

void MeshElementRemoval::searchByMaterialID(const std::vector<MeshLib::Element*> & ele_vec, unsigned matID)
{
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getValue()==matID)
			matchedIDs.push_back(i);
		i++;
	}
	this->updateUnion(matchedIDs);
}

void MeshElementRemoval::searchByElementType(const std::vector<MeshLib::Element*> & ele_vec, MeshElemType eleType)
{
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getGeomType()==eleType)
			matchedIDs.push_back(i);
		i++;
	}
	this->updateUnion(matchedIDs);
}

void MeshElementRemoval::searchByZeroContent(const std::vector<MeshLib::Element*> & ele_vec)
{
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getContent()<std::numeric_limits<double>::epsilon())
			matchedIDs.push_back(i);
		i++;
	}
	this->updateUnion(matchedIDs);
}

void MeshElementRemoval::searchByBoundingBox(const std::vector<MeshLib::Element*> & ele_vec, const MeshLib::Node &x1, const MeshLib::Node &x2)
{
	std::vector<MeshLib::Node> extent;
	extent.push_back(x1); extent.push_back(x2);
	const GeoLib::AABB<MeshLib::Node> aabb(extent.begin(), extent.end());

	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) 
	{
		std::size_t nElemNodes (ele_vec[i]->getNNodes());
		for (std::size_t j=0; j<nElemNodes; ++j)
			if (!aabb.containsPoint(*ele_vec[i]->getNode(j)))
			{
				matchedIDs.push_back(i);
				break;
			}
		i++;
	}
	this->updateUnion(matchedIDs);
}

void MeshElementRemoval::updateUnion(const std::vector<std::size_t> &vec)
{
	std::vector<std::size_t> vec_temp(vec.size() + _marked_elements.size());
	auto it = std::set_union(vec.begin(), vec.end(), _marked_elements.begin(), _marked_elements.end(), vec_temp.begin());
	vec_temp.resize(it - vec_temp.begin());
	_marked_elements.assign(vec_temp.begin(), vec_temp.end());
}

std::vector<MeshLib::Element*> MeshElementRemoval::excludeElements(const std::vector<MeshLib::Element*> & vec_src_eles, const std::vector<std::size_t> &vec_removed) const
{
	std::vector<MeshLib::Element*> vec_dest_eles(vec_src_eles.size() - vec_removed.size());
	std::size_t k=0;
	for (std::size_t i=0; i<vec_src_eles.size(); i++) {
		if (std::find(vec_removed.begin(), vec_removed.end(), i) == vec_removed.end()) {
			vec_dest_eles[k] = vec_src_eles[i];
			k++;
		}
	}
	return vec_dest_eles;
}

void MeshElementRemoval::copyNodesElements(	const std::vector<MeshLib::Node*> &src_nodes, const std::vector<MeshLib::Element*> &src_elems,
						std::vector<MeshLib::Node*> &dst_nodes, std::vector<MeshLib::Element*> &dst_elems) const 
{
	// copy nodes
	dst_nodes.resize(src_nodes.size());
	for (std::size_t i=0; i<dst_nodes.size(); i++) {
		dst_nodes[i] = new MeshLib::Node(*src_nodes[i]);
	}

	// copy elements with new nodes
	dst_elems.resize(src_elems.size());
	for (std::size_t i=0; i<dst_elems.size(); i++) {
		auto* src_elem = src_elems[i];
		auto* dst_elem = src_elem->clone();
		for (unsigned j=0; j<src_elem->getNNodes(); j++) {
			dst_elem->setNode(j, dst_nodes[src_elem->getNode(j)->getID()]);
		}
		dst_elems[i] = dst_elem;
	}
}





} // end namespace MeshLib
