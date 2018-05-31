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

/// A subset of nodes or elements on a single mesh.
class MeshSubset
{
public:
    /// Construct a mesh subset from vector of nodes on the given mesh.
    /// \param msh Mesh
    /// \param vec_items Vector of Node pointers.
    /// \param delete_ptr Deletes the vector of Node pointers if true.
    /// \note When delete_ptr is set only the vector is deleted, not the
    /// elements of the vector.
    MeshSubset(const Mesh& msh, std::vector<Node*> const* vec_items,
        bool const delete_ptr = false)
        : _msh(msh), _nodes(vec_items), _eles(nullptr), _delete_ptr(delete_ptr)
    {}

    /// Construct a mesh subset from vector of elements on the given mesh.
    /// \param msh Mesh
    /// \param vec_items Vector of Element pointers.
    /// \param delete_ptr Deletes the vector of Element pointers if true.
    /// \note When delete_ptr is set only the vector is deleted, not the
    /// elements of the vector.
    MeshSubset(const Mesh& msh, std::vector<Element*> const* vec_items,
        bool const delete_ptr = false)
        : _msh(msh), _nodes(nullptr), _eles(vec_items), _delete_ptr(delete_ptr)
    {}

    /// construct from both nodes and elements
    /// Construct a mesh subset from vector of nodes and a vector of elements on
    /// the given mesh.
    /// \param msh Mesh
    /// \param vec_nodes Vector of Node pointers.
    /// \param vec_eles Vector of Element pointers.
    /// \param delete_ptr Deletes the vector of Node pointers if true.
    /// \note When delete_ptr is set only the vectors are deleted, not the
    /// elements of the vectors.
    MeshSubset(const Mesh& msh, std::vector<Node*> const* vec_nodes,
        std::vector<Element*> const* vec_eles, bool const delete_ptr = false)
        : _msh(msh), _nodes(vec_nodes), _eles(vec_eles), _delete_ptr(delete_ptr)
    {}

    ~MeshSubset()
    {
        if (_delete_ptr)
        {
            delete _nodes;
            delete _eles;
        }
    }

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

    /// Constructs a new mesh subset which is a set intersection of the current
    /// nodes and the provided vector of nodes.
    /// An empty mesh subset may be returned, not a nullptr, in case of empty
    /// intersection or empty input vector.
    MeshSubset getIntersectionByNodes(std::vector<Node*> const& nodes) const
    {
        auto* active_nodes = new std::vector<Node*>;

        if (_nodes == nullptr || _nodes->empty())
            return MeshSubset(_msh, active_nodes);  // Empty mesh subset

        for (auto n : nodes)
        {
            auto it = std::find(_nodes->cbegin(), _nodes->cend(), n);
            if (it == _nodes->cend())
                continue;
            active_nodes->push_back(n);
        }

        // Transfer the ownership of active_nodes to the new MeshSubset, which
        // deletes the pointer itself.
        return MeshSubset(_msh, active_nodes,
                          false);  // This causes a memory leak of the non
                                   // deleted active_nodes vector.
                                   // Calling ctor with 'true', causes double
                                   // free of the nodes vector, when MS is
                                   // copied multiple times.
    }

    Mesh const& getMesh() const
    {
        return _msh;
    }

private:
    Mesh const& _msh;
    std::vector<Node*> const* _nodes;
    std::vector<Element*> const* _eles;
    bool const _delete_ptr = false;

};

}   // namespace MeshLib
