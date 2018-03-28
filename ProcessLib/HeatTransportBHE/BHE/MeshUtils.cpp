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
void getBHEDataInMesh(
    MeshLib::Mesh const& mesh,
    std::vector<MeshLib::Element*>& vec_soil_elements,
    std::vector<int>& vec_BHE_mat_IDs,
    std::vector<std::vector<MeshLib::Element*>>& vec_BHE_elements,
    std::vector<MeshLib::Node*>& vec_soil_nodes,
    std::vector<std::vector<MeshLib::Node*>>& vec_BHE_nodes)
{
    // get vectors of matrix elements and BHE elements
    vec_soil_elements.reserve(mesh.getNumberOfElements());
    std::vector<MeshLib::Element*> all_BHE_elements;
    for (MeshLib::Element* e : mesh.getElements())
    {
        // As for the first step, all elements are counted as a soil element
        // first. Those elements connected with a BHE will picked up and
        // reorganized into a seperate vector at the end of the function.
        if (e->getDimension() == mesh.getDimension())
            vec_soil_elements.push_back(e);
        else if (e->getDimension() == (unsigned int)1)
            all_BHE_elements.push_back(e);
    }

    // get BHE material IDs
    auto opt_material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));
    for (MeshLib::Element* e : all_BHE_elements)
        vec_BHE_mat_IDs.push_back((*opt_material_ids)[e->getID()]);
    BaseLib::makeVectorUnique(vec_BHE_mat_IDs);
    DBUG("-> found %d BHE material groups", vec_BHE_mat_IDs.size());

    // create a vector of BHE elements for each group
    vec_BHE_elements.resize(vec_BHE_mat_IDs.size());
    for (unsigned bhe_id = 0; bhe_id < vec_BHE_mat_IDs.size(); bhe_id++)
    {
        const auto bhe_mat_id = vec_BHE_mat_IDs[bhe_id];
        std::vector<MeshLib::Element*>& vec_elements = vec_BHE_elements[bhe_id];
        std::copy_if(all_BHE_elements.begin(), all_BHE_elements.end(),
                     std::back_inserter(vec_elements),
                     [&](MeshLib::Element* e) {
                         return (*opt_material_ids)[e->getID()] == bhe_mat_id;
                     });
        DBUG("-> found %d elements on the BHE_%d", vec_elements.size(), bhe_id);
    }

    // get a vector of BHE nodes
    vec_BHE_nodes.resize(vec_BHE_mat_IDs.size());
    for (unsigned bhe_id = 0; bhe_id < vec_BHE_mat_IDs.size(); bhe_id++)
    {
        std::vector<MeshLib::Node*>& vec_nodes = vec_BHE_nodes[bhe_id];
        for (MeshLib::Element* e : vec_BHE_elements[bhe_id])
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
    vec_soil_nodes.reserve(mesh.getNumberOfNodes());
    for (MeshLib::Node* n : mesh.getNodes())
    {
        // All elements are counted as a soil element
        vec_soil_nodes.push_back(n);
    }

    // final count of 3 types of elements
    // They are
    // (i)  soil,
    // (ii) BHE
    DBUG("-> found total %d soil elements and %d BHE elements",
         vec_soil_elements.size(),
         all_BHE_elements.size());
}
}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
