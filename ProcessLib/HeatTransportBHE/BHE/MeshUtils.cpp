/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshUtils.h"

#include "BaseLib/Algorithm.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
BHEMeshData getBHEDataInMesh(MeshLib::Mesh const& mesh)
{
    BHEMeshData bheMeshData;
    // partition all the mesh elements, and copy them into
    // two seperate vectors, one with only matrix elements
    // and the other only BHE elements
    bheMeshData.soil_elements.reserve(mesh.getNumberOfElements());
    std::vector<MeshLib::Element*> all_BHE_elements;
    auto& all_mesh_elements = mesh.getElements();
    std::partition_copy(
        std::begin(all_mesh_elements),
        std::end(all_mesh_elements),
        std::back_inserter(all_BHE_elements),
        std::back_inserter(bheMeshData.soil_elements),
        [](MeshLib::Element* e) { return e->getDimension() == 1; });

    // get BHE material IDs
    auto opt_material_ids = MeshLib::materialIDs(mesh);
    if (opt_material_ids == nullptr)
    {
        OGS_FATAL("Not able to get material IDs! ");
    }
    for (MeshLib::Element* e : all_BHE_elements)
        bheMeshData.BHE_mat_IDs.push_back((*opt_material_ids)[e->getID()]);
    BaseLib::makeVectorUnique(bheMeshData.BHE_mat_IDs);
    DBUG("-> found %d BHE material groups", bheMeshData.BHE_mat_IDs.size());

    // create a vector of BHE elements for each group
    bheMeshData.BHE_elements.resize(bheMeshData.BHE_mat_IDs.size());
    for (unsigned bhe_id = 0; bhe_id < bheMeshData.BHE_mat_IDs.size(); bhe_id++)
    {
        const auto bhe_mat_id = bheMeshData.BHE_mat_IDs[bhe_id];
        std::vector<MeshLib::Element*>& vec_elements =
            bheMeshData.BHE_elements[bhe_id];
        std::copy_if(all_BHE_elements.begin(), all_BHE_elements.end(),
                     std::back_inserter(vec_elements),
                     [&](MeshLib::Element* e) {
                         return (*opt_material_ids)[e->getID()] == bhe_mat_id;
                     });
        DBUG("-> found %d elements on the BHE_%d", vec_elements.size(), bhe_id);
    }

    // get a vector of BHE nodes
    bheMeshData.BHE_nodes.resize(bheMeshData.BHE_mat_IDs.size());
    for (unsigned bhe_id = 0; bhe_id < bheMeshData.BHE_mat_IDs.size(); bhe_id++)
    {
        std::vector<MeshLib::Node*>& vec_nodes = bheMeshData.BHE_nodes[bhe_id];
        for (MeshLib::Element* e : bheMeshData.BHE_elements[bhe_id])
        {
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                vec_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
            }
        }
        BaseLib::makeVectorUnique(
            vec_nodes, [](MeshLib::Node* node1, MeshLib::Node* node2) {
                return node1->getID() < node2->getID();
            });
        DBUG("-> found %d nodes on the BHE_%d", vec_nodes.size(), bhe_id);
    }

    // get all the nodes into the pure soil nodes vector
    bheMeshData.soil_nodes.reserve(mesh.getNumberOfNodes());
    for (MeshLib::Node* n : mesh.getNodes())
    {
        // All elements are counted as a soil element
        bheMeshData.soil_nodes.push_back(n);
    }

    // finalLy counting two types of elements
    // They are (i) soil, and (ii) BHE type of elements
    DBUG("-> found total %d soil elements and %d BHE elements",
         bheMeshData.soil_elements.size(),
         bheMeshData.BHE_elements.size());

    return bheMeshData;
}
}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib