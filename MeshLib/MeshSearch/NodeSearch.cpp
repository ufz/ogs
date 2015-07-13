/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NodeSearch.h"

#include <set>

#include <logog/include/logog.hpp>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib {

NodeSearch::NodeSearch(const MeshLib::Mesh &mesh)
	: _mesh(mesh)
{
}

NodeSearch::~NodeSearch()
{
}

std::size_t NodeSearch::searchByElementIDs(const std::vector<std::size_t> &elements, bool only_match_all_connected_elements)
{
	std::vector<std::size_t> connected_nodes;
	if (only_match_all_connected_elements)
	{
		std::vector<std::size_t> node_marked_counts(_mesh.getNNodes(), 0); //this approach is not optimum for memory size
		std::for_each(elements.begin(), elements.end(),
			[&](std::size_t eid)
			{
				auto* e = _mesh.getElement(eid);
				for (unsigned i=0; i<e->getNBaseNodes(); i++) {
					node_marked_counts[e->getNodeIndex(i)]++;
				}
			});
		for (std::size_t i=0; i<node_marked_counts.size(); i++)
		{
			if (node_marked_counts[i] == _mesh.getNode(i)->getElements().size())
				connected_nodes.push_back(i);
		}
	} else {
		std::for_each(elements.begin(), elements.end(),
			[&](std::size_t eid)
			{
				auto* e = _mesh.getElement(eid);
				for (unsigned i=0; i<e->getNBaseNodes(); i++) {
					connected_nodes.push_back(e->getNodeIndex(i));
				}
			});
		std::sort(connected_nodes.begin(), connected_nodes.end());
		auto it = std::unique(connected_nodes.begin(), connected_nodes.end());
		connected_nodes.resize(std::distance(connected_nodes.begin(),it));
	}
	this->updateUnion(connected_nodes);
	return connected_nodes.size();
}

void NodeSearch::updateUnion(const std::vector<std::size_t> &vec)
{
	std::vector<std::size_t> vec_temp(vec.size() + _marked_nodes.size());
	auto it = std::set_union(vec.begin(), vec.end(), _marked_nodes.begin(), _marked_nodes.end(), vec_temp.begin());
	vec_temp.resize(it - vec_temp.begin());
	_marked_nodes.assign(vec_temp.begin(), vec_temp.end());
}

std::vector<Node*> getNodes(std::vector<Element*> const& elements)
{
	std::set<Node*> nodes_set;
	for (auto e : elements)
	{
		Node* const* nodes = e->getNodes();
		unsigned const nnodes = e->getNNodes();
		nodes_set.insert(nodes, nodes + nnodes);
	}

	std::vector<Node*> nodes;
	nodes.reserve(nodes_set.size());

	std::move(nodes_set.cbegin(), nodes_set.cend(),
		std::back_inserter(nodes));

	return nodes;
}

} // end namespace MeshLib
