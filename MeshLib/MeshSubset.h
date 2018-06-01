/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <vector>

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"


namespace MeshLib
{
inline std::vector<Node*> nodesNodesIntersection(
    std::vector<Node*> const& nodes_a, std::vector<Node*> const& nodes_b)
{
    if (nodes_a.empty() || nodes_b.empty())
    {
        return {};
    }

    std::vector<Node*> active_nodes;

    for (auto const& n_a : nodes_a)
    {
        auto it = std::find(begin(nodes_b), end(nodes_b), n_a);
        if (it != end(nodes_b))
        {
            active_nodes.push_back(n_a);
        }
    }

    return active_nodes;
}

/// A subset of nodes or elements on a single mesh.
class MeshSubset
{
public:
    /// Construct a mesh subset from vector of nodes on the given mesh.
    /// \param msh Mesh
    /// \param vec_items Vector of Node pointers.
    MeshSubset(const Mesh& msh, std::vector<Node*> const* vec_items)
        : _msh(msh), _nodes(vec_items), _eles(nullptr)
    {}

    /// Construct a mesh subset from vector of elements on the given mesh.
    /// \param msh Mesh
    /// \param vec_items Vector of Element pointers.
    MeshSubset(const Mesh& msh, std::vector<Element*> const* vec_items)
        : _msh(msh), _nodes(nullptr), _eles(vec_items)
    {}

    /// construct from both nodes and elements
    /// Construct a mesh subset from vector of nodes and a vector of elements on
    /// the given mesh.
    /// \param msh Mesh
    /// \param vec_nodes Vector of Node pointers.
    /// \param vec_eles Vector of Element pointers.
    MeshSubset(const Mesh& msh, std::vector<Node*> const* vec_nodes,
               std::vector<Element*> const* vec_eles)
        : _msh(msh), _nodes(vec_nodes), _eles(vec_eles)
    {}

    /// return the total number of mesh items
    std::size_t getNumberOfTotalItems() const
    {
        return getNumberOfNodes() + getNumberOfElements();
    }

    /// return this mesh ID
    std::size_t getMeshID() const
    {
        return _msh.getID();
    }

    /// return the number of registered nodes
    std::size_t getNumberOfNodes() const
    {
        return (_nodes==nullptr) ? 0 : _nodes->size();
    }

    /// Returns the global node id Node::getID() of i-th node in the mesh
    /// subset.
    /// \pre The _nodes must be a valid pointer to a vector of size > i.
    std::size_t getNodeID(std::size_t const i) const
    {
        assert(_nodes && i < _nodes->size());
        return (*_nodes)[i]->getID();
    }

    /// return the number of registered elements
    std::size_t getNumberOfElements() const
    {
        return (_eles==nullptr) ? 0 : _eles->size();
    }

    /// Returns the global element id Element::getID() of i-th element in the
    /// mesh subset.
    /// \pre The _eles must be a valid pointer to a vector of size > i.
    std::size_t getElementID(std::size_t const i) const
    {
        assert(_eles && i < _eles->size());
        return (*_eles)[i]->getID();
    }

    std::vector<Element*>::const_iterator elementsBegin() const
    {
        return _msh.getElements().cbegin();
    }

    std::vector<Element*>::const_iterator elementsEnd() const
    {
        return _msh.getElements().cend();
    }

    std::vector<Node*> const& getNodes() const
    {
        assert(_nodes);
        return *_nodes;
    }

    Mesh const& getMesh() const
    {
        return _msh;
    }

private:
    Mesh const& _msh;
    std::vector<Node*> const* _nodes;
    std::vector<Element*> const* _eles;
};
}   // namespace MeshLib
