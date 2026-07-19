// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshUtils.h"

#include <algorithm>
#include <set>
#include <unordered_map>
#include <utility>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"

namespace
{
std::vector<MeshLib::Element*> extractOneDimensionalElements(
    std::vector<MeshLib::Element*> const& elements)
{
    std::vector<MeshLib::Element*> one_dimensional_elements;

    copy_if(
        begin(elements), end(elements), back_inserter(one_dimensional_elements),
        [](MeshLib::Element const* const e) { return e->getDimension() == 1; });

    return one_dimensional_elements;
}

std::vector<int> getUniqueMaterialIds(
    MeshLib::PropertyVector<int> const& material_ids,
    std::vector<MeshLib::Element*> const& elements)
{
    std::set<int> unique_material_ids;
    std::transform(begin(elements), end(elements),
                   inserter(unique_material_ids, end(unique_material_ids)),
                   [&material_ids](MeshLib::Element const* const e)
                   { return material_ids[e->getID()]; });
    return {begin(unique_material_ids), end(unique_material_ids)};
}
}  // namespace

namespace ProcessLib
{
namespace HeatTransportBHE
{
BHEMeshData getBHEDataInMesh(MeshLib::Mesh const& mesh)
{
    std::vector<MeshLib::Element*> const all_bhe_elements =
        extractOneDimensionalElements(mesh.getElements());

    // finally counting two types of elements
    // They are (i) soil, and (ii) BHE type of elements
    DBUG("-> found total {:d} soil elements and {:d} BHE elements",
         mesh.getNumberOfElements() - all_bhe_elements.size(),
         all_bhe_elements.size());

    // get BHE material IDs
    auto const* const opt_material_ids = MeshLib::materialIDs(mesh);
    if (opt_material_ids == nullptr)
    {
        OGS_FATAL("Not able to get material IDs! ");
    }
    auto const& material_ids = *opt_material_ids;

    auto const& bhe_material_ids =
        getUniqueMaterialIds(material_ids, all_bhe_elements);
    DBUG("-> found {:d} BHE material groups", bhe_material_ids.size());

    // create a vector of BHE elements for each group
    std::vector<std::vector<MeshLib::Element*>> bhe_elements;
    bhe_elements.resize(bhe_material_ids.size());
    for (unsigned bhe_id = 0; bhe_id < bhe_material_ids.size(); bhe_id++)
    {
        const auto bhe_mat_id = bhe_material_ids[bhe_id];
        std::vector<MeshLib::Element*>& vec_elements = bhe_elements[bhe_id];
        copy_if(begin(all_bhe_elements), end(all_bhe_elements),
                back_inserter(vec_elements),
                [&](MeshLib::Element const* const e)
                { return material_ids[e->getID()] == bhe_mat_id; });
        DBUG("-> found {:d} elements on the BHE_{:d}", vec_elements.size(),
             bhe_id);
    }

    // get a vector of BHE nodes
    std::vector<std::vector<MeshLib::Node*>> bhe_nodes;
    bhe_nodes.resize(bhe_material_ids.size());
    for (unsigned bhe_id = 0; bhe_id < bhe_material_ids.size(); bhe_id++)
    {
        std::vector<MeshLib::Node*>& vec_nodes = bhe_nodes[bhe_id];
        for (MeshLib::Element* e : bhe_elements[bhe_id])
        {
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                vec_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
            }
        }
        BaseLib::makeVectorUnique(vec_nodes,
                                  MeshLib::idsComparator<MeshLib::Node*>);

        DBUG("-> found {:d} nodes on the BHE_{:d}", vec_nodes.size(), bhe_id);
    }

    return {bhe_material_ids, bhe_elements, bhe_nodes};
}

BHEEndpoints findBHEEndpointsFromElementOrdering(
    std::vector<MeshLib::Element*> const& bhe_elements)
{
    if (bhe_elements.empty())
    {
        OGS_FATAL(
            "findBHEEndpointsFromElementOrdering called with an empty BHE "
            "element list.");
    }

    // Adjacency: node id -> list of (element index, local-node-index 0 or 1)
    // that reference it. For a well-formed BHE chain every node appears in
    // either 1 (endpoint) or 2 (interior) such records.
    std::unordered_map<std::size_t, std::vector<std::pair<std::size_t, int>>>
        adjacency;
    for (std::size_t ei = 0; ei < bhe_elements.size(); ++ei)
    {
        auto const* e = bhe_elements[ei];
        if (e->getNumberOfNodes() != 2)
        {
            OGS_FATAL(
                "BHE element {:d} has {:d} nodes; expected a 2-node line "
                "element.",
                e->getID(), e->getNumberOfNodes());
        }
        for (int local = 0; local < 2; ++local)
        {
            adjacency[e->getNode(local)->getID()].emplace_back(ei, local);
        }
    }

    // Endpoints are nodes connected to exactly one BHE line element. Any
    // node connected to more than 2 elements means the BHE is branching and
    // is not a simple chain.
    std::vector<std::size_t> endpoint_node_ids;
    for (auto const& [nid, refs] : adjacency)
    {
        if (refs.size() == 1)
        {
            endpoint_node_ids.push_back(nid);
        }
        else if (refs.size() > 2)
        {
            OGS_FATAL(
                "BHE element chain has a branching junction at node {:d} "
                "(connected to {:d} elements). Expected a simple line chain.",
                nid, refs.size());
        }
    }

    if (endpoint_node_ids.size() != 2)
    {
        OGS_FATAL(
            "BHE elements do not form a single connected line chain (found "
            "{:d} endpoints; expected 2). Check that all line elements with "
            "the same BHE material ID share endpoints node-to-node.",
            endpoint_node_ids.size());
    }

    // The inlet is the endpoint that appears as local node 0 of its element;
    // the outlet is the endpoint that appears as local node 1 of its
    // element. If neither endpoint is at local node 0 (or both are), the
    // chain-start or chain-end element is reversed.
    MeshLib::Node const* inlet = nullptr;
    MeshLib::Node const* outlet = nullptr;
    for (std::size_t const nid : endpoint_node_ids)
    {
        auto const& [ei, local] = adjacency[nid].front();
        if (local == 0)
        {
            if (inlet != nullptr)
            {
                OGS_FATAL(
                    "Both BHE chain endpoints (nodes {:d} and {:d}) appear "
                    "as local node 0 of their connected element; chain "
                    "orientation is ambiguous. Reorder element nodes so the "
                    "chain runs node 0 -> node 1 from inlet to outlet.",
                    inlet->getID(), nid);
            }
            inlet = bhe_elements[ei]->getNode(0);
        }
        else
        {
            outlet = bhe_elements[ei]->getNode(1);
        }
    }

    if (inlet == nullptr || outlet == nullptr)
    {
        OGS_FATAL(
            "BHE chain endpoints are not consistent with the elements' "
            "node-0 -> node-1 ordering. Reorder element nodes so that the "
            "intended inlet endpoint is node 0 of the chain-start element "
            "and the intended outlet endpoint is node 1 of the chain-end "
            "element.");
    }

    // Walk the chain from the inlet: at each step, the unused element whose
    // node 0 equals the current node must exist, otherwise some interior
    // element is reversed.
    std::vector<bool> used(bhe_elements.size(), false);
    MeshLib::Node const* cur = inlet;
    for (std::size_t step = 0; step < bhe_elements.size(); ++step)
    {
        auto const& refs = adjacency[cur->getID()];
        bool advanced = false;
        for (auto const& [ei, local] : refs)
        {
            if (!used[ei] && local == 0)
            {
                used[ei] = true;
                cur = bhe_elements[ei]->getNode(1);
                advanced = true;
                break;
            }
        }
        if (!advanced)
        {
            OGS_FATAL(
                "BHE element chain is inconsistently oriented at node {:d}: "
                "an interior element does not have this node as its node 0. "
                "Reorder the element's nodes so the chain runs node 0 -> "
                "node 1 from inlet to outlet.",
                cur->getID());
        }
    }

    if (cur != outlet)
    {
        OGS_FATAL(
            "BHE chain walk did not terminate at the expected outlet "
            "endpoint (expected node {:d}, reached node {:d}).",
            outlet->getID(), cur->getID());
    }

    return {inlet, outlet};
}

}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
