/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NODESEARCH_H_
#define NODESEARCH_H_

#include <vector>

namespace MeshLib
{

// forward declarations
class Mesh;
class Element;
class Node;

/// Node search class
class NodeSearch final
{
public:
	explicit NodeSearch(const MeshLib::Mesh &mesh);

	/// return marked node IDs
	const std::vector<std::size_t>& getSearchedNodeIDs() const {return _marked_nodes; }

	/// Marks all nodes connecting to any of the given elements
	std::size_t searchByElementIDs(const std::vector<std::size_t> &element_ids, bool only_match_all_connected_elements = false)
	{
		std::vector<std::size_t> connected_nodes =
			(only_match_all_connected_elements
				? searchByElementIDsMatchAllConnectedElements(element_ids)
				: searchByElementIDsNotMatchAllConnectedElements(element_ids));

		this->updateUnion(connected_nodes);
		return connected_nodes.size();
	}

	/// Marks all unused nodes
	std::size_t searchUnused();

private:
	std::vector<std::size_t> searchByElementIDsMatchAllConnectedElements(const std::vector<std::size_t> &element_ids);
	std::vector<std::size_t> searchByElementIDsNotMatchAllConnectedElements(const std::vector<std::size_t> &element_ids);

	/// Updates the vector of marked items with values from vec.
	void updateUnion(const std::vector<std::size_t> &vec);

	/// The mesh from which elements should be removed.
	const MeshLib::Mesh &_mesh;
	/// The vector of element indices that should be removed.
	std::vector<std::size_t> _marked_nodes;
};

/// Create a vector of unique nodes used by given elements.
std::vector<Node*> getUniqueNodes(std::vector<Element*> const& elements);

} // end namespace MeshLib

#endif //NODESEARCH_H_
