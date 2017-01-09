
/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "AppendLinesAlongPolyline.h"

#include <logog/include/logog.hpp>

#include "GeoLib/Polyline.h"
#include "GeoLib/PolylineVec.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"

#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"

namespace MeshGeoToolsLib
{
std::unique_ptr<MeshLib::Mesh> appendLinesAlongPolylines(
    const MeshLib::Mesh& mesh, const GeoLib::PolylineVec& ply_vec)
{
    // copy existing nodes and elements
    std::vector<MeshLib::Node*> vec_new_nodes = MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> vec_new_eles = MeshLib::copyElementVector(mesh.getElements(), vec_new_nodes);

    std::vector<int> new_mat_ids;
    {
        auto mat_ids = mesh.getProperties().getPropertyVector<int>("MaterialIDs");
        if (mat_ids) {
            new_mat_ids.reserve((*mat_ids).size());
            std::copy((*mat_ids).cbegin(), (*mat_ids).cend(),
                std::back_inserter(new_mat_ids));
        }
    }
    int max_matID(0);
    if (!new_mat_ids.empty())
        max_matID = *(std::max_element(new_mat_ids.cbegin(), new_mat_ids.cend()));

    const std::size_t n_ply (ply_vec.size());
    // for each polyline
    for (std::size_t k(0); k < n_ply; k++)
    {
        const GeoLib::Polyline* ply = (*ply_vec.getVector())[k];

        // search nodes on the polyline
        MeshGeoToolsLib::MeshNodesAlongPolyline mshNodesAlongPoly(mesh, *ply, mesh.getMinEdgeLength()*0.5);
        auto &vec_nodes_on_ply = mshNodesAlongPoly.getNodeIDs();
        if (vec_nodes_on_ply.empty()) {
            std::string ply_name;
            ply_vec.getNameOfElementByID(k, ply_name);
            INFO("No nodes found on polyline %s", ply_name.c_str());
            continue;
        }

        // add line elements
        for (std::size_t i=0; i<vec_nodes_on_ply.size()-1; i++) {
            std::array<MeshLib::Node*, 2> element_nodes;
            element_nodes[0] = vec_new_nodes[vec_nodes_on_ply[i]];
            element_nodes[1] = vec_new_nodes[vec_nodes_on_ply[i+1]];
            vec_new_eles.push_back(
                new MeshLib::Line(element_nodes, vec_new_eles.size()));
            new_mat_ids.push_back(max_matID+k+1);
        }
    }

    // generate a mesh
    const std::string name = mesh.getName() + "_with_lines";
    std::unique_ptr<MeshLib::Mesh> new_mesh(
        new MeshLib::Mesh(name, vec_new_nodes, vec_new_eles));
    auto opt_mat_pv = new_mesh->getProperties().createNewPropertyVector<int>(
        "MaterialIDs", MeshLib::MeshItemType::Cell);
    if (opt_mat_pv) {
        auto & mat_pv = *opt_mat_pv;
        mat_pv.reserve(new_mat_ids.size());
        std::copy(new_mat_ids.cbegin(), new_mat_ids.cend(),
            std::back_inserter(mat_pv));
    }
    return new_mesh;
}

} // MeshGeoToolsLib

