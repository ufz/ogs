/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
    /// \param use_taylor_hood_elements Flag to indicate whether the Taylor-Hood
    ///                                 elements are used.
    MeshSubset(const Mesh& msh, std::vector<Node*> const& vec_items,
               const bool use_taylor_hood_elements = false)
        : _msh(msh),
          _nodes(vec_items),
          _use_taylor_hood_elements(use_taylor_hood_elements)
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
            auto node_is_part_of_mesh =
                [&mesh_nodes = sorted_nodes](MeshLib::Node* const& n)
            {
                auto it = lower_bound(begin(mesh_nodes), end(mesh_nodes), n);
                if (it == end(mesh_nodes))
                {
                    ERR("A node {:d} ({:g}, {:g}, {:g}) in mesh subset is not "
                        "a part of the mesh.",
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
    std::size_t getMeshID() const { return _msh.getID(); }

    bool useTaylorHoodElements() const { return _use_taylor_hood_elements; }

    std::vector<Element*>::const_iterator elementsBegin() const
    {
        return _msh.getElements().cbegin();
    }

    std::vector<Element*>::const_iterator elementsEnd() const
    {
        return _msh.getElements().cend();
    }

    std::vector<Node*> const& getNodes() const { return _nodes; }

    Mesh const& getMesh() const { return _msh; }

private:
    Mesh const& _msh;
    std::vector<Node*> const& _nodes;
    bool const _use_taylor_hood_elements;
};

namespace views
{
inline auto meshLocations(MeshSubset const& mesh_subset,
                          MeshItemType const item_type)
{
    return mesh_subset.getNodes() | ids |
           ranges::views::transform(
               [mesh_id = mesh_subset.getMeshID(),
                item_type](std::size_t const node_id)
               { return Location{mesh_id, item_type, node_id}; });
}
}  // namespace views
}  // namespace MeshLib
