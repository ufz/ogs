/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHSUBSET_H_
#define MESHSUBSET_H_

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
    /// Throws std::out_of_range exception if there are no nodes available.
    std::size_t getNodeID(std::size_t const i) const
    {
        if (!_nodes)
            throw std::out_of_range(
                "In MeshSubset::getNodeID(): no nodes or nodes are empty.");

        return (*_nodes)[i]->getID();
    }

    /// return the number of registered elements
    std::size_t getNElements() const
    {
        return (_eles==nullptr) ? 0 : _eles->size();
    }

    std::vector<Element*>::const_iterator elementsBegin() const
    {
        return _msh.getElements().cbegin();
    }

    std::vector<Element*>::const_iterator elementsEnd() const
    {
        return _msh.getElements().cend();
    }

private:
    const Mesh& _msh;
    std::vector<Node*> const* _nodes;
    std::vector<Element*> const* _eles;

};

}   // namespace MeshLib

#endif  // MESHSUBSET_H_
