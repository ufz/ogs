/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "AppendLinesAlongPolyline.h"

#include <range/v3/view/any_view.hpp>
#include <range/v3/view/common.hpp>
#include <range/v3/view/repeat_n.hpp>

#include "BaseLib/Logging.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/PolylineVec.h"
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"

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
    assert(material_ids != nullptr);
    int const max_matID = material_ids ? ranges::max(*material_ids) : 0;

    auto const edgeLengths = minMaxEdgeLength(mesh.getElements());
    double const min_edge = edgeLengths.first;

    std::vector<int> new_mat_ids;
    const std::size_t n_ply(ply_vec.size());
    // for each polyline
    for (std::size_t k(0); k < n_ply; k++)
    {
        auto const* const ply = ply_vec.getVector()[k];

        // search nodes on the polyline
        MeshGeoToolsLib::MeshNodesAlongPolyline mshNodesAlongPoly(
            mesh, *ply, min_edge * 0.5, MeshGeoToolsLib::SearchAllNodes::Yes);
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
        std::make_unique<MeshLib::Mesh>(name, vec_new_nodes, vec_new_eles,
                                        true /* compute_element_neighbors */);
    auto new_material_ids =
        new_mesh->getProperties().createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell);
    if (!new_material_ids)
    {
        OGS_FATAL("Could not create MaterialIDs cell vector in new mesh.");
    }

    auto initial_values =
        material_ids ? ranges::any_view<int const>(*material_ids)
                     : ranges::views::repeat_n(0, mesh.getNumberOfElements());

    new_material_ids->assign(ranges::views::common(
        ranges::views::concat(initial_values, new_mat_ids)));

    return new_mesh;
}

}  // namespace MeshGeoToolsLib
