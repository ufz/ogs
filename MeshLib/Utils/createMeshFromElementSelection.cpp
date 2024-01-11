/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "createMeshFromElementSelection.h"

#include <range/v3/numeric.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/map.hpp>
#include <unordered_map>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/addPropertyToMesh.h"

namespace MeshLib
{
std::unique_ptr<MeshLib::Mesh> createMeshFromElementSelection(
    std::string mesh_name, std::vector<MeshLib::Element*> const& elements)
{
    auto ids_vector = views::ids | ranges::to<std::vector>();

    DBUG("Found {:d} elements in the mesh", elements.size());

    // Store bulk element ids for each of the new elements.
    auto bulk_element_ids = elements | ids_vector;

    // original node ids to newly created nodes.
    std::unordered_map<std::size_t, MeshLib::Node*> id_node_hash_map;
    id_node_hash_map.reserve(
        elements.size());  // There will be at least one node per element.

    for (auto& e : elements)
    {
        // For each node find a cloned node in map or create if there is none.
        unsigned const n_nodes = e->getNumberOfNodes();
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            const MeshLib::Node* n = e->getNode(i);
            auto const it = id_node_hash_map.find(n->getID());
            if (it == id_node_hash_map.end())
            {
                auto new_node_in_map = id_node_hash_map[n->getID()] =
                    new MeshLib::Node(*n);
                e->setNode(i, new_node_in_map);
            }
            else
            {
                e->setNode(i, it->second);
            }
        }
    }

    std::map<std::size_t, MeshLib::Node*> nodes_map;
    for (const auto& n : id_node_hash_map)
    {
        nodes_map[n.first] = n.second;
    }

    // Copy the unique nodes pointers.
    auto element_nodes =
        nodes_map | ranges::views::values | ranges::to<std::vector>;

    // Store bulk node ids for each of the new nodes.
    auto bulk_node_ids =
        nodes_map | ranges::views::keys | ranges::to<std::vector>;

    auto mesh = std::make_unique<MeshLib::Mesh>(
        std::move(mesh_name), std::move(element_nodes), std::move(elements),
        true /* compute_element_neighbors */);
    assert(mesh != nullptr);

    addPropertyToMesh(*mesh, getBulkIDString(MeshLib::MeshItemType::Cell),
                      MeshLib::MeshItemType::Cell, 1, bulk_element_ids);
    addPropertyToMesh(*mesh, getBulkIDString(MeshLib::MeshItemType::Node),
                      MeshLib::MeshItemType::Node, 1, bulk_node_ids);

    return mesh;
}

}  // namespace MeshLib
