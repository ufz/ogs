// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshUtils.h"

#include <algorithm>
#include <array>
#include <set>
#include <unordered_map>

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

std::array<MeshLib::Node*, 2> getElementEndpoints(MeshLib::Element& element)
{
    if (element.getGeomType() != MeshLib::MeshElemType::LINE)
    {
        OGS_FATAL(
            "Expected a line element for BHE element id {:d}, but got {:s}.",
            element.getID(),
            MeshLib::MeshElemType2String(element.getGeomType()));
    }

    return {element.getNode(0), element.getNode(1)};
}

/// Build a node-to-elements adjacency map for a set of BHE elements.
std::unordered_map<std::size_t, std::vector<MeshLib::Element*>>
buildNodeToElementMap(std::vector<MeshLib::Element*> const& elements)
{
    std::unordered_map<std::size_t, std::vector<MeshLib::Element*>> map;
    for (auto* element : elements)
    {
        auto const [n0, n1] = getElementEndpoints(*element);
        map[n0->getID()].push_back(element);
        map[n1->getID()].push_back(element);
    }
    return map;
}

/// Find the wellhead node by filtering BHE nodes.
/// Record the nodes that are connected to only one element.
/// Then verifying that there are exactly two nodes exist,
/// and returning the one with the greater z-coordinate.
MeshLib::Node* findWellheadNode(
    std::vector<MeshLib::Node*> const& bhe_nodes,
    std::unordered_map<std::size_t, std::vector<MeshLib::Element*>> const&
        node_to_elements)
{
    if (bhe_nodes.empty())
    {
        OGS_FATAL("BHE node list is empty while constructing mesh data.");
    }

    std::vector<MeshLib::Node*> endpoints;
    std::copy_if(
        bhe_nodes.begin(), bhe_nodes.end(), std::back_inserter(endpoints),
        [&node_to_elements](MeshLib::Node* node)
        {
            auto const it = node_to_elements.find(node->getID());
            return it != node_to_elements.end() && it->second.size() == 1;
        });

    if (endpoints.size() != 2)
    {
        OGS_FATAL(
            "The BHE mesh must form a single continuous linear chain with "
            "exactly 2 endpoints. Found {:d} endpoints.",
            endpoints.size());
    }

    auto* a = endpoints[0];
    auto* b = endpoints[1];

    if ((*a)[2] == (*b)[2])
    {
        OGS_FATAL(
            "Both BHE chain endpoints share the same z-coordinate ({:g}). "
            "The wellhead choice is ambiguous.",
            (*a)[2]);
    }

    return ((*a)[2] > (*b)[2]) ? a : b;
}

/// Walk the BHE element chain starting from the wellhead node using mesh
/// topology (node-element adjacency).  Returns elements in order from wellhead
/// to bottom and computes arc-length element distances.
///
/// Unlike z-sorting, this works correctly for any BHE orientation (vertical,
/// inclined, or horizontal).
struct ChainWalkResult
{
    std::vector<MeshLib::Element*> ordered_elements;
    std::vector<MeshLib::Node*> ordered_nodes;
    std::unordered_map<std::size_t, double> element_distances_from_wellhead;
};

ChainWalkResult walkChainFromWellhead(
    std::vector<MeshLib::Element*> const& bhe_elements,
    std::vector<MeshLib::Node*> const& bhe_nodes)
{
    if (bhe_elements.empty())
    {
        return {};
    }

    auto const node_to_elements = buildNodeToElementMap(bhe_elements);
    MeshLib::Node* wellhead = findWellheadNode(bhe_nodes, node_to_elements);

    // The wellhead node must be an endpoint of exactly one element.
    auto const it = node_to_elements.find(wellhead->getID());
    if (it == node_to_elements.end() || it->second.empty())
    {
        OGS_FATAL("Wellhead node {:d} is not connected to any BHE element.",
                  wellhead->getID());
    }
    if (it->second.size() != 1)
    {
        OGS_FATAL(
            "Wellhead node {:d} is connected to {:d} BHE elements; "
            "expected exactly 1 (chain endpoint).",
            wellhead->getID(), it->second.size());
    }

    std::vector<MeshLib::Element*> ordered;
    ordered.reserve(bhe_elements.size());
    std::vector<MeshLib::Node*> ordered_nodes;
    ordered_nodes.reserve(bhe_elements.size() + 1);
    ordered_nodes.push_back(wellhead);
    std::unordered_map<std::size_t, double> distances;
    distances.reserve(bhe_elements.size());

    MeshLib::Node* prev_node = wellhead;
    MeshLib::Element* current = it->second.front();
    double accumulated_distance = 0.0;

    while (current != nullptr)
    {
        ordered.push_back(current);
        double const len = current->computeVolume();
        distances[current->getID()] = accumulated_distance + 0.5 * len;
        accumulated_distance += len;

        // Find exit node (the endpoint not shared with prev_node).
        auto const [n0, n1] = getElementEndpoints(*current);
        MeshLib::Node* exit_node = (n0 == prev_node) ? n1 : n0;
        ordered_nodes.push_back(exit_node);

        // Find the next element connected to exit_node.
        auto const exit_it = node_to_elements.find(exit_node->getID());
        MeshLib::Element* next = nullptr;
        if (exit_it != node_to_elements.end())
        {
            auto const& elems = exit_it->second;
            auto const next_it =
                std::find_if(elems.begin(), elems.end(),
                             [current](auto* e) { return e != current; });
            if (next_it != elems.end())
            {
                next = *next_it;
            }
        }

        prev_node = exit_node;
        current = next;
    }

    if (ordered.size() != bhe_elements.size())
    {
        OGS_FATAL(
            "BHE chain walk visited {:d} elements but the group has {:d}. "
            "The BHE mesh must form a single continuous linear chain.",
            ordered.size(), bhe_elements.size());
    }

    return {std::move(ordered), std::move(ordered_nodes), std::move(distances)};
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

    std::unordered_map<std::size_t, double> bhe_element_distances_from_wellhead;
    std::vector<std::vector<MeshLib::Node*>> bhe_topology_ordered_nodes(
        bhe_material_ids.size());

    std::unordered_map<std::size_t, int> bhe_element_section_indices;

    for (unsigned bhe_id = 0; bhe_id < bhe_material_ids.size(); bhe_id++)
    {
        auto walk_result =
            walkChainFromWellhead(bhe_elements[bhe_id], bhe_nodes[bhe_id]);

        bhe_elements[bhe_id] = std::move(walk_result.ordered_elements);
        bhe_topology_ordered_nodes[bhe_id] =
            std::move(walk_result.ordered_nodes);

        for (auto const& [element_id, distance] :
             walk_result.element_distances_from_wellhead)
        {
            bhe_element_distances_from_wellhead[element_id] = distance;
        }
    }

    return {bhe_material_ids,
            bhe_elements,
            bhe_nodes,
            bhe_topology_ordered_nodes,
            std::move(bhe_element_distances_from_wellhead),
            std::move(bhe_element_section_indices)};
}

void BHEMeshData::updateElementSectionIndices(
    std::vector<BHE::BHETypes> const& bhes)
{
    if (bhes.size() != BHE_elements.size())
    {
        OGS_FATAL(
            "Mismatch between number of BHE definitions ({:d}) and BHE mesh "
            "groups ({:d}) while creating per-element section indices.",
            bhes.size(), BHE_elements.size());
    }

    for (std::size_t bhe_id = 0; bhe_id < bhes.size(); ++bhe_id)
    {
        for (auto const* const element : BHE_elements[bhe_id])
        {
            auto const element_id = element->getID();
            auto const dist_it =
                BHE_element_distances_from_wellhead.find(element_id);
            if (dist_it == BHE_element_distances_from_wellhead.end())
            {
                OGS_FATAL(
                    "Could not determine section index: no distance found for "
                    "element id {:d}.",
                    element_id);
            }

            double const distance_from_wellhead = dist_it->second;

            BHE_element_section_indices[element_id] = visit(
                [distance_from_wellhead](auto const& bhe)
                {
                    return bhe.borehole_geometry.sections.getSectionIndex(
                        distance_from_wellhead);
                },
                bhes[bhe_id]);
        }
    }
}

}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
