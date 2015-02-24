/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHSUBSET_H_
#define MESHSUBSET_H_

#include <cassert>
#include <vector>

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"


namespace MeshLib
{

/// A subset of nodes or elements on a single mesh.
class MeshSubset
{
public:
    /// construct from nodes
    MeshSubset(const Mesh& msh, std::vector<Node*> const& vec_items)
        : _msh(msh), _nodes(&vec_items), _eles(nullptr)
    {}

    /// construct from elements
    MeshSubset(const Mesh& msh, std::vector<Element*> const& vec_items)
        : _msh(msh), _nodes(nullptr), _eles(&vec_items)
    {}

    /// construct from both nodes and elements
    MeshSubset(const Mesh& msh, std::vector<Node*> const& vec_nodes,
               std::vector<Element*> const& vec_eles)
        : _msh(msh), _nodes(&vec_nodes), _eles(&vec_eles)
    {}

    ~MeshSubset() {}

    /// return the total number of mesh items
    std::size_t getNTotalItems() const
    {
        return getNNodes() + getNElements();
    }

    /// return this mesh ID
    std::size_t getMeshID() const
    {
        return _msh.getID();
    }

    /// return the number of registered nodes
    std::size_t getNNodes() const
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
    std::size_t getNElements() const
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

    /// Constructs a new mesh subset which is a set intersection of the current
    /// nodes and the provided vector of nodes.
    /// An empty mesh subset may be returned, not a nullptr, in case of empty
    /// intersection or empty input vector.
    MeshSubset*
    getIntersectionByNodes(std::vector<Node*> const& nodes) const
    {
        std::vector<Node*>* active_nodes = new std::vector<Node*>;

        if (_nodes == nullptr || _nodes->empty())
            return new MeshSubset(_msh, *active_nodes);   // Empty mesh subset

        for (auto n : nodes)
        {
            auto it = std::find(_nodes->cbegin(), _nodes->cend(), n);
            if (it == _nodes->cend())
                continue;
            active_nodes->push_back(n);
        }

        return new MeshSubset(_msh, *active_nodes);
    }

private:
    const Mesh& _msh;
    std::vector<Node*> const* _nodes;
    std::vector<Element*> const* _eles;

};

}   // namespace MeshLib

#endif  // MESHSUBSET_H_
