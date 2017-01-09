/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

    /// Marks all nodes connected to any of the given elements ids.
    /// \return number of connected nodes.
    std::size_t searchNodesConnectedToOnlyGivenElements(const std::vector<std::size_t> &element_ids);

    /// Marks all unused nodes
    std::size_t searchUnused();

    /// Marks all boundary nodes
    std::size_t searchBoundaryNodes();

private:
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
