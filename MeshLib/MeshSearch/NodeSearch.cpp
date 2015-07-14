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

std::vector<std::size_t> NodeSearch::searchByElementIDsMatchAllConnectedElements(const std::vector<std::size_t> &elements)
{
	std::vector<std::size_t> connected_nodes;

	std::vector<std::size_t> node_marked_counts(_mesh.getNNodes(), 0);
	//this approach is not optimal for memory size

	for(std::size_t eid : elements)
	{
		auto* e = _mesh.getElement(eid);
		for (unsigned i=0; i<e->getNBaseNodes(); i++) {
			node_marked_counts[e->getNodeIndex(i)]++;
		}
	}

	for (std::size_t i=0; i<node_marked_counts.size(); i++)
	{
		if (node_marked_counts[i] == _mesh.getNode(i)->getElements().size())
			connected_nodes.push_back(i);
	}
	return connected_nodes;
}

std::vector<std::size_t> NodeSearch::searchByElementIDsNotMatchAllConnectedElements(const std::vector<std::size_t> &elements)
{
	std::vector<std::size_t> connected_nodes;

	for(std::size_t eid : elements)
	{
		auto* e = _mesh.getElement(eid);
		for (unsigned i=0; i<e->getNBaseNodes(); i++) {
			connected_nodes.push_back(e->getNodeIndex(i));
		}
	}

	std::sort(connected_nodes.begin(), connected_nodes.end());
	auto it = std::unique(connected_nodes.begin(), connected_nodes.end());
	connected_nodes.resize(std::distance(connected_nodes.begin(),it));

	return connected_nodes;
}

std::size_t NodeSearch::searchUnused()
{
	const size_t nNodes (_mesh.getNNodes());
	const std::vector<MeshLib::Node*> &nodes (_mesh.getNodes());
	std::vector<std::size_t> del_node_idx;

	for (unsigned i=0; i<nNodes; ++i)
		if (nodes[i]->getNElements() == 0)
			del_node_idx.push_back(i);

	this->updateUnion(del_node_idx);
	return del_node_idx.size();
}

void NodeSearch::updateUnion(const std::vector<std::size_t> &vec)
{
	std::vector<std::size_t> vec_temp(vec.size() + _marked_nodes.size());
	auto it = std::set_union(vec.begin(), vec.end(), _marked_nodes.begin(), _marked_nodes.end(), vec_temp.begin());
	vec_temp.resize(it - vec_temp.begin());
	_marked_nodes.assign(vec_temp.begin(), vec_temp.end());
}

std::vector<Node*> getUniqueNodes(std::vector<Element*> const& elements)
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
