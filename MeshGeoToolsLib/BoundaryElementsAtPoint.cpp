/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BoundaryElementsAtPoint.h"

#include "GeoLib/Point.h"
#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Point.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/Node.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{
BoundaryElementsAtPoint::BoundaryElementsAtPoint(
    MeshLib::Mesh const& mesh, MeshNodeSearcher const& mshNodeSearcher,
    GeoLib::Point const& point, const bool multiple_nodes_allowed)
    : mesh_(mesh), point_(point)
{
    auto const node_ids = mshNodeSearcher.getMeshNodeIDs(point_);
    if (node_ids.empty())
    {
        OGS_FATAL(
            "BoundaryElementsAtPoint: the mesh node searcher was unable to "
            "locate the point ({:f}, {:f}, {:f}) in the mesh.",
            point_[0], point_[1], point_[2]);
    }
    if (node_ids.size() == 1)
    {
        std::array<MeshLib::Node*, 1> const nodes = {
            {const_cast<MeshLib::Node*>(mesh_.getNode(node_ids[0]))}};
        boundary_elements_.push_back(new MeshLib::Point{nodes, node_ids[0]});
        return;
    }

    auto& mesh_nodes =
        const_cast<std::vector<MeshLib::Node*>&>(mesh.getNodes());
    std::size_t const nearest_node_id =
        *std::min_element(node_ids.begin(), node_ids.end(),
                          [&mesh_nodes, &point](auto id0, auto id1) {
                              return MathLib::sqrDist(*mesh_nodes[id0], point) <
                                     MathLib::sqrDist(*mesh_nodes[id1], point);
                          });

    if (!multiple_nodes_allowed)
    {
        OGS_FATAL(
            "BoundaryElementsAtPoint: the mesh node searcher found {:d} points "
            "near the requested point ({:f}, {:f}, {:f}) in the mesh, while "
            "exactly one is expected. Node  (id={:d}) ({:f}, {:f}, {:f}) has "
            "distance {:f}.",
            node_ids.size(), point_[0], point_[1], point_[2],
            mesh_nodes[nearest_node_id]->getID(),
            (*mesh_nodes[nearest_node_id])[0],
            (*mesh_nodes[nearest_node_id])[1],
            (*mesh_nodes[nearest_node_id])[2],
            MathLib::sqrDist(*mesh_nodes[nearest_node_id], point));
    }
    WARN(
        "BoundaryElementsAtPoint: the mesh node searcher found {:d} points "
        "near the requested point ({:f}, {:f}, {:f}) in the mesh, while "
        "exactly one is expected. Node  (id={:d}) ({:f}, {:f}, {:f}) has "
        "distance {:f}.",
        node_ids.size(), point_[0], point_[1], point_[2],
        mesh_nodes[nearest_node_id]->getID(), (*mesh_nodes[nearest_node_id])[0],
        (*mesh_nodes[nearest_node_id])[1], (*mesh_nodes[nearest_node_id])[2],
        MathLib::sqrDist(*mesh_nodes[nearest_node_id], point));

    std::array<MeshLib::Node*, 1> const nodes = {
        {const_cast<MeshLib::Node*>(mesh_.getNode(nearest_node_id))}};

    boundary_elements_.push_back(
        new MeshLib::Point{nodes, nearest_node_id});
}

BoundaryElementsAtPoint::~BoundaryElementsAtPoint()
{
    for (auto p : boundary_elements_)
    {
        delete p;
    }
}
}  // namespace MeshGeoToolsLib
