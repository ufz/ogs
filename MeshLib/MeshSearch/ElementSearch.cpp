/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementSearch.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib {

ElementSearch::ElementSearch(const MeshLib::Mesh &mesh)
	: _mesh(mesh)
{
}

ElementSearch::~ElementSearch()
{
}

std::size_t ElementSearch::searchByMaterialID(int const matID)
{
	boost::optional<MeshLib::PropertyVector<int> const&> opt_pv(
		this->_mesh.getProperties().getPropertyVector<int>("MaterialIDs")
	);
	if (!opt_pv)
		return 0;

	MeshLib::PropertyVector<int> const& pv(opt_pv.get());

	std::vector<std::size_t> matchedIDs;
	for (std::size_t i(0); i<pv.getNumberOfTuples(); ++i) {
		if (pv[i]==matID)
			matchedIDs.push_back(i);
	}

	this->updateUnion(matchedIDs);
	return matchedIDs.size();
}

std::size_t ElementSearch::searchByElementType(MeshElemType eleType)
{
	const std::vector<MeshLib::Element*> &ele_vec (this->_mesh.getElements());
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getGeomType()==eleType)
			matchedIDs.push_back(i);
		i++;
	}
	this->updateUnion(matchedIDs);
	return matchedIDs.size();
}

std::size_t ElementSearch::searchByContent(double eps)
{
	const std::vector<MeshLib::Element*> &ele_vec (this->_mesh.getElements());
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getContent() < eps)
			matchedIDs.push_back(i);
		i++;
	}
	this->updateUnion(matchedIDs);
	return matchedIDs.size();
}

std::size_t ElementSearch::searchByBoundingBox(
	GeoLib::AABB<MathLib::Point3d> const& aabb)
{
	const std::vector<MeshLib::Element*> &ele_vec (this->_mesh.getElements());

	std::vector<std::size_t> matchedIDs;
	const std::size_t n_elems(ele_vec.size());
	for (std::size_t i = 0; i<n_elems; i++)
	{
		std::size_t nElemNodes (ele_vec[i]->getNBaseNodes());
		for (std::size_t j=0; j<nElemNodes; ++j)
			if (!aabb.containsPoint(*ele_vec[i]->getNode(j)))
			{
				matchedIDs.push_back(i);
				break;
			}
	}
	this->updateUnion(matchedIDs);
	return matchedIDs.size();
}

std::size_t ElementSearch::searchByNodeIDs(const std::vector<std::size_t> &nodes)
{
	std::vector<std::size_t> connected_elements;
	std::for_each(nodes.begin(), nodes.end(),
		[&](std::size_t node_id)
		{
			for (auto* e : _mesh.getNode(node_id)->getElements()) {
				connected_elements.push_back(e->getID());
			}
		});
	std::sort(connected_elements.begin(), connected_elements.end());
	auto it = std::unique(connected_elements.begin(), connected_elements.end());
	connected_elements.resize(std::distance(connected_elements.begin(),it));
	this->updateUnion(connected_elements);
	return connected_elements.size();
}

void ElementSearch::updateUnion(const std::vector<std::size_t> &vec)
{
	std::vector<std::size_t> vec_temp(vec.size() + _marked_elements.size());
	auto it = std::set_union(vec.begin(), vec.end(), _marked_elements.begin(), _marked_elements.end(), vec_temp.begin());
	vec_temp.resize(it - vec_temp.begin());
	_marked_elements.assign(vec_temp.begin(), vec_temp.end());
}

} // end namespace MeshLib
