// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
    std::size_t searchNodesConnectedToOnlyGivenElements(
        const std::vector<std::size_t>& elements);

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
