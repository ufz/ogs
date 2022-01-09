/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshUtils.h"

#include <set>

#include "BaseLib/Algorithm.h"
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

    copy_if(begin(elements), end(elements),
            back_inserter(one_dimensional_elements),
            [](MeshLib::Element* e) { return e->getDimension() == 1; });

    return one_dimensional_elements;
}

std::vector<int> getUniqueMaterialIds(
    std::vector<int> const& material_ids,
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
                [&](MeshLib::Element* e)
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
                                  [](MeshLib::Node* node1, MeshLib::Node* node2)
                                  { return node1->getID() < node2->getID(); });
        DBUG("-> found {:d} nodes on the BHE_{:d}", vec_nodes.size(), bhe_id);
    }

    return {bhe_material_ids, bhe_elements, bhe_nodes};
}
}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
