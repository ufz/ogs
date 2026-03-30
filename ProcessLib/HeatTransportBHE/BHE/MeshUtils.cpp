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

std::array<MeshLib::Node const*, 2> getElementEndpoints(
    MeshLib::Element const& element)
{
    if (element.getNumberOfBaseNodes() < 2)
    {
        OGS_FATAL(
            "Expected a line element with at least two base nodes for BHE "
            "element id {:d}, but got {:d} base nodes.",
            element.getID(), element.getNumberOfBaseNodes());
    }

    return {element.getNode(0), element.getNode(1)};
}

/// Sorts BHE elements from wellhead (highest z) to bottom (lowest z) so that
/// a simple linear chain walk gives arc-length distances in order.
void sortBheElementsFromWellhead(std::vector<MeshLib::Element*>& elements)
{
    std::sort(
        elements.begin(), elements.end(),
        [](MeshLib::Element const* const a, MeshLib::Element const* const b)
        {
            auto const za = ((*a->getNode(0))[2] + (*a->getNode(1))[2]) / 2.0;
            auto const zb = ((*b->getNode(0))[2] + (*b->getNode(1))[2]) / 2.0;
            return za > zb;
        });
}

MeshLib::Node const& findWellheadNode(
    std::vector<MeshLib::Element*> const& sorted_elements)
{
    if (sorted_elements.empty())
    {
        OGS_FATAL("BHE element list is empty while constructing mesh data.");
    }

    // After sorting, the wellhead is the higher-z endpoint of the first
    // element.
    auto const [n0, n1] = getElementEndpoints(*sorted_elements.front());
    return (*n0)[2] >= (*n1)[2] ? *n0 : *n1;
}

std::unordered_map<std::size_t, double> computeElementDistancesFromWellhead(
    std::vector<MeshLib::Element*> const& sorted_bhe_elements)
{
    std::unordered_map<std::size_t, double> distances;
    distances.reserve(sorted_bhe_elements.size());

    if (sorted_bhe_elements.empty())
    {
        return distances;
    }

    // BHE elements form a linear chain. Walking from the wellhead element by
    // element gives arc-length distances in a single O(N) pass.
    MeshLib::Node const* prev_node = &findWellheadNode(sorted_bhe_elements);
    double accumulated_distance = 0.0;

    for (auto* const element : sorted_bhe_elements)
    {
        auto const [n0, n1] = getElementEndpoints(*element);
        // The exit node is the endpoint not shared with the previous element.
        MeshLib::Node const* const bottom = (n0 == prev_node) ? n1 : n0;
        double const len = element->computeVolume();
        distances[element->getID()] = accumulated_distance + 0.5 * len;
        accumulated_distance += len;
        prev_node = bottom;
    }

    return distances;
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

    std::unordered_map<std::size_t, int> bhe_element_section_indices;

    for (unsigned bhe_id = 0; bhe_id < bhe_material_ids.size(); bhe_id++)
    {
        sortBheElementsFromWellhead(bhe_elements[bhe_id]);

        for (auto const& [element_id, distance] :
             computeElementDistancesFromWellhead(bhe_elements[bhe_id]))
        {
            bhe_element_distances_from_wellhead[element_id] = distance;
        }
    }

    return {bhe_material_ids, bhe_elements, bhe_nodes,
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
