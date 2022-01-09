/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "AppendLinesAlongPolyline.h"

#include "BaseLib/Logging.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/PolylineVec.h"
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
std::unique_ptr<MeshLib::Mesh> appendLinesAlongPolylines(
    const MeshLib::Mesh& mesh, const GeoLib::PolylineVec& ply_vec)
{
    // copy existing nodes and elements
    std::vector<MeshLib::Node*> vec_new_nodes =
        MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> vec_new_eles =
        MeshLib::copyElementVector(mesh.getElements(), vec_new_nodes);

    auto const material_ids = materialIDs(mesh);
    int const max_matID =
        material_ids
            ? *(std::max_element(begin(*material_ids), end(*material_ids)))
            : 0;

    std::vector<int> new_mat_ids;
    const std::size_t n_ply(ply_vec.size());
    // for each polyline
    for (std::size_t k(0); k < n_ply; k++)
    {
        const GeoLib::Polyline* ply = (*ply_vec.getVector())[k];

        // search nodes on the polyline
        MeshGeoToolsLib::MeshNodesAlongPolyline mshNodesAlongPoly(
            mesh, *ply, mesh.getMinEdgeLength() * 0.5,
            MeshGeoToolsLib::SearchAllNodes::Yes);
        auto& vec_nodes_on_ply = mshNodesAlongPoly.getNodeIDs();
        if (vec_nodes_on_ply.empty())
        {
            std::string ply_name;
            ply_vec.getNameOfElementByID(k, ply_name);
            INFO("No nodes found on polyline {:s}", ply_name);
            continue;
        }

        // add line elements
        for (std::size_t i = 0; i < vec_nodes_on_ply.size() - 1; i++)
        {
            std::array<MeshLib::Node*, 2> element_nodes;
            element_nodes[0] = vec_new_nodes[vec_nodes_on_ply[i]];
            element_nodes[1] = vec_new_nodes[vec_nodes_on_ply[i + 1]];
            vec_new_eles.push_back(
                new MeshLib::Line(element_nodes, vec_new_eles.size()));
            new_mat_ids.push_back(max_matID + k + 1);
        }
    }

    // generate a mesh
    const std::string name = mesh.getName() + "_with_lines";
    auto new_mesh =
        std::make_unique<MeshLib::Mesh>(name, vec_new_nodes, vec_new_eles);
    auto new_material_ids =
        new_mesh->getProperties().createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell);
    if (!new_material_ids)
    {
        OGS_FATAL("Could not create MaterialIDs cell vector in new mesh.");
    }
    new_material_ids->reserve(new_mesh->getNumberOfElements());
    if (material_ids != nullptr)
    {
        std::copy(begin(*material_ids), end(*material_ids),
                  std::back_inserter(*new_material_ids));
    }
    else
    {
        new_material_ids->resize(mesh.getNumberOfElements());
    }
    std::copy(begin(new_mat_ids), end(new_mat_ids),
              std::back_inserter(*new_material_ids));
    return new_mesh;
}

}  // namespace MeshGeoToolsLib
