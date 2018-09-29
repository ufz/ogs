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
/// A subset of nodes on a single mesh.
class MeshSubset
{
public:
    /// Construct a mesh subset from vector of nodes on the given mesh.
    /// \param msh Mesh
    /// \param vec_items Vector of Node pointers.
    MeshSubset(const Mesh& msh, std::vector<Node*> const& vec_items)
        : _msh(msh), _nodes(vec_items)
    {
        // If the mesh nodes and the given nodes point to the same vector, they
        // must be equal.
        if (&_msh.getNodes() == &_nodes)
        {
            return;
        }

        //
        // Testing if the given nodes belong to the mesh.
        //
        {
            // Need sorted version of the large vector.
            auto sorted_nodes = _msh.getNodes();  // full copy of pointers.
            sort(begin(sorted_nodes), end(sorted_nodes));

            // Then proceed with the search function.
            auto node_is_part_of_mesh = [& mesh_nodes = sorted_nodes](
                                            MeshLib::Node* const& n) {
                auto it = lower_bound(begin(mesh_nodes), end(mesh_nodes), n);
                if (it == end(mesh_nodes))
                {
                    ERR("A node %d (%g, %g, %g) in mesh subset is not a part "
                        "of the mesh.",
                        n->getID(), (*n)[0], (*n)[1], (*n)[2]);
                    return false;
                }
                return true;
            };
            if (!std::all_of(begin(_nodes), end(_nodes), node_is_part_of_mesh))
            {
                OGS_FATAL("The mesh subset construction failed.");
            }
        }
    }

    /// return this mesh ID
    std::size_t getMeshID() const
    {
        return _msh.getID();
    }

    /// return the number of registered nodes
    std::size_t getNumberOfNodes() const
    {
        return _nodes.size();
    }

    /// Returns the global node id Node::getID() of i-th node in the mesh
    /// subset.
    /// \pre The _nodes vector must be of size > i.
    std::size_t getNodeID(std::size_t const i) const
    {
        assert(i < _nodes.size());
        return _nodes[i]->getID();
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
        return _nodes;
    }

    Mesh const& getMesh() const
    {
        return _msh;
    }

private:
    Mesh const& _msh;
    std::vector<Node*> const& _nodes;
};
}   // namespace MeshLib
