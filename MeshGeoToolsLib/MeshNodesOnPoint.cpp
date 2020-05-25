/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshNodesOnPoint.h"

#include "MeshLib/Mesh.h"

namespace MeshGeoToolsLib
{
MeshNodesOnPoint::MeshNodesOnPoint(MeshLib::Mesh const& mesh,
                                   GeoLib::Grid<MeshLib::Node> const& mesh_grid,
                                   GeoLib::Point const& pnt,
                                   double epsilon_radius,
                                   SearchAllNodes search_all_nodes)
    : mesh_(mesh), pnt_(pnt)
{
    std::vector<std::size_t> vec_ids(
        mesh_grid.getPointsInEpsilonEnvironment(pnt, epsilon_radius));
    if (search_all_nodes == SearchAllNodes::Yes)
    {
        msh_node_ids_ = vec_ids;
    }
    else
    {
        for (auto id : vec_ids)
        {
            if (mesh.isBaseNode(id))
            {
                msh_node_ids_.push_back(id);
            }
        }
    }
}

}  // end namespace MeshGeoToolsLib

