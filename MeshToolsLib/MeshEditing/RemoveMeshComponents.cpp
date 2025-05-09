/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RemoveMeshComponents.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"

namespace MeshToolsLib
{
MeshLib::Mesh* removeElements(
    const MeshLib::Mesh& mesh,
    const std::vector<std::size_t>& removed_element_ids,
    const std::string& new_mesh_name)
{
    if (removed_element_ids.empty())
    {
        INFO("No elements to remove");
        return nullptr;
    }

    INFO("Removing total {:d} elements...", removed_element_ids.size());
    std::vector<MeshLib::Element*> tmp_elems =
        BaseLib::excludeObjectCopy(mesh.getElements(), removed_element_ids);
    INFO("{:d} elements remain in mesh.", tmp_elems.size());

    // copy node and element objects
    std::vector<MeshLib::Node*> new_nodes =
        MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> new_elems =
        MeshLib::copyElementVector(tmp_elems, new_nodes);

    // delete unused nodes
    MeshLib::NodeSearch ns(mesh);
    ns.searchNodesConnectedToOnlyGivenElements(removed_element_ids);
    auto& removed_node_ids(ns.getSearchedNodeIDs());
    INFO("Removing total {:d} nodes...", removed_node_ids.size());
    for (auto nodeid : removed_node_ids)
    {
        delete new_nodes[nodeid];
        new_nodes[nodeid] = nullptr;
    }
    new_nodes.erase(std::remove(new_nodes.begin(), new_nodes.end(), nullptr),
                    new_nodes.end());

    if (!new_elems.empty())
    {
        MeshLib::Mesh* new_mesh =
            new MeshLib::Mesh(new_mesh_name, new_nodes, new_elems,
                              true /* compute_element_neighbors */,
                              mesh.getProperties().excludeCopyProperties(
                                  removed_element_ids, removed_node_ids));
        return new_mesh;
    }

    INFO("Current selection removes all elements.");
    return nullptr;
}

std::vector<bool> markUnusedNodes(
    std::vector<MeshLib::Element*> const& elements,
    std::vector<MeshLib::Node*> const& nodes)
{
    std::vector<bool> unused_nodes(nodes.size(), true);
    for (auto e : elements)
    {
        for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
        {
            unused_nodes[getNodeIndex(*e, i)] = false;
        }
    }

    return unused_nodes;
}

void removeMarkedNodes(std::vector<bool> const& nodes_to_delete,
                       std::vector<MeshLib::Node*>& nodes)
{
    assert(nodes_to_delete.size() == nodes.size());

    for (std::size_t i = 0; i < nodes.size(); i++)
    {
        if (nodes_to_delete[i])
        {
            delete nodes[i];
            nodes[i] = nullptr;
        }
    }
    nodes.erase(remove(begin(nodes), end(nodes), nullptr), end(nodes));
}

MeshLib::Mesh* removeNodes(const MeshLib::Mesh& mesh,
                           const std::vector<std::size_t>& del_nodes_idx,
                           const std::string& new_mesh_name)
{
    if (del_nodes_idx.empty())
    {
        return nullptr;
    }

    // copy node and element objects
    std::vector<MeshLib::Node*> new_nodes =
        MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> new_elems =
        MeshLib::copyElementVector(mesh.getElements(), new_nodes);

    // delete elements
    MeshLib::ElementSearch es(mesh);
    es.searchByNodeIDs(del_nodes_idx);
    auto& removed_element_ids = es.getSearchedElementIDs();
    for (auto eid : removed_element_ids)
    {
        delete new_elems[eid];
        new_elems[eid] = nullptr;
    }
    new_elems.erase(std::remove(new_elems.begin(), new_elems.end(), nullptr),
                    new_elems.end());

    // check unused nodes due to element deletion
    std::vector<bool> const node_delete_flag =
        markUnusedNodes(new_elems, new_nodes);

    // delete unused nodes
    removeMarkedNodes(node_delete_flag, new_nodes);

    if (!new_elems.empty())
    {
        MeshLib::Mesh* new_mesh =
            new MeshLib::Mesh(new_mesh_name, new_nodes, new_elems,
                              true /* compute_element_neighbors */,
                              mesh.getProperties().excludeCopyProperties(
                                  removed_element_ids, del_nodes_idx));
        return new_mesh;
    }

    return nullptr;
}
}  // namespace MeshToolsLib
