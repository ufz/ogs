/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RemoveMeshComponents.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "DuplicateMeshComponents.h"

namespace MeshLib
{

namespace details
{

std::vector<MeshLib::Element*> excludeElementCopy(
    std::vector<MeshLib::Element*> const& vec_src_eles,
    std::vector<std::size_t> const& vec_removed)
{
    std::vector<MeshLib::Element*> vec_dest_eles(vec_src_eles.size()-vec_removed.size());

    unsigned cnt (0);
    for (std::size_t i=0; i<vec_removed[0]; ++i)
        vec_dest_eles[cnt++] = vec_src_eles[i];
    for (std::size_t i=1; i<vec_removed.size(); ++i)
        for (std::size_t j=vec_removed[i-1]+1; j<vec_removed[i]; ++j)
            vec_dest_eles[cnt++] = vec_src_eles[j];
    for (std::size_t i=vec_removed.back()+1; i<vec_src_eles.size(); ++i)
        vec_dest_eles[cnt++] = vec_src_eles[i];

    return vec_dest_eles;
}

} // details

MeshLib::Mesh* removeElements(const MeshLib::Mesh& mesh, const std::vector<std::size_t> &removed_element_ids, const std::string &new_mesh_name)
{
    if (removed_element_ids.empty())
    {
        INFO("No elements to remove");
        return nullptr;
    }

    INFO("Removing total %d elements...", removed_element_ids.size());
    std::vector<MeshLib::Element*> tmp_elems = details::excludeElementCopy(
        mesh.getElements(),
        removed_element_ids
    );
    INFO("%d elements remain in mesh.", tmp_elems.size());

    // copy node and element objects
    std::vector<MeshLib::Node*> new_nodes = MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> new_elems = MeshLib::copyElementVector(tmp_elems, new_nodes);

    // delete unused nodes
    NodeSearch ns(mesh);
    ns.searchNodesConnectedToOnlyGivenElements(removed_element_ids);
    auto &removed_node_ids(ns.getSearchedNodeIDs());
    INFO("Removing total %d nodes...", removed_node_ids.size());
    for (auto nodeid : removed_node_ids)
    {
        delete new_nodes[nodeid];
        new_nodes[nodeid] = nullptr;
    }
    new_nodes.erase(std::remove(new_nodes.begin(), new_nodes.end(), nullptr), new_nodes.end());

    if (!new_elems.empty())
    {
        MeshLib::Mesh* new_mesh =
            new MeshLib::Mesh(new_mesh_name, new_nodes, new_elems,
                              mesh.getProperties().excludeCopyProperties(
                                  removed_element_ids, removed_node_ids));
        return new_mesh;
    }
    else
    {
        INFO("Current selection removes all elements.");
        return nullptr;
    }
}

MeshLib::Mesh* removeNodes(const MeshLib::Mesh &mesh, const std::vector<std::size_t> &del_nodes_idx, const std::string &new_mesh_name)
{
    if (del_nodes_idx.empty())
        return nullptr;

    // copy node and element objects
    std::vector<MeshLib::Node*> new_nodes = MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> new_elems = MeshLib::copyElementVector(mesh.getElements(), new_nodes);

    // delete elements
    MeshLib::ElementSearch es(mesh);
    es.searchByNodeIDs(del_nodes_idx);
    auto& removed_element_ids = es.getSearchedElementIDs();
    for (auto eid : removed_element_ids)
    {
        delete new_elems[eid];
        new_elems[eid] = nullptr;
    }
    new_elems.erase(std::remove(new_elems.begin(), new_elems.end(), nullptr), new_elems.end());

    // check unused nodes due to element deletion
    std::vector<bool> node_delete_flag(new_nodes.size(), true);
    for (auto e : new_elems) {
        for (unsigned i=0; i<e->getNumberOfNodes(); i++)
            node_delete_flag[e->getNodeIndex(i)] = false;
    }

    // delete unused nodes
    for (std::size_t i=0; i<new_nodes.size(); i++)
    {
        if (!node_delete_flag[i]) continue;
        delete new_nodes[i];
        new_nodes[i] = nullptr;
    }
    new_nodes.erase(std::remove(new_nodes.begin(), new_nodes.end(), nullptr), new_nodes.end());

    if (!new_elems.empty())
    {
        MeshLib::Mesh* new_mesh =
            new MeshLib::Mesh(new_mesh_name, new_nodes, new_elems,
                              mesh.getProperties().excludeCopyProperties(
                                  removed_element_ids, del_nodes_idx));
        return new_mesh;
    }
    else
    {
        return nullptr;
    }
}
} // end namespace MeshLib

