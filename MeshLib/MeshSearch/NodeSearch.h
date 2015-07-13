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
class NodeSearch
{
public:
	explicit NodeSearch(const MeshLib::Mesh &mesh);

	~NodeSearch();

	/// return marked node IDs
	const std::vector<std::size_t>& getSearchedNodeIDs() const {return _marked_nodes; }

	/// Marks all nodes connecting to any of the given elements
	std::size_t searchByElementIDs(const std::vector<std::size_t> &element_ids, bool only_match_all_connected_elements = false);

private:
	/// Updates the vector of marked items with values from vec.
	void updateUnion(const std::vector<std::size_t> &vec);

	/// The mesh from which elements should be removed.
	const MeshLib::Mesh &_mesh;
	/// The vector of element indices that should be removed.
	std::vector<std::size_t> _marked_nodes;
};

/// Create a vector of unique nodes used by given elements.
std::vector<Node*> getNodes(std::vector<Element*> const& elements);

} // end namespace MeshLib

#endif //NODESEARCH_H_
