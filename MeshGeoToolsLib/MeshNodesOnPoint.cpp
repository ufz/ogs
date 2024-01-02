/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshNodesOnPoint.h"

#include <range/v3/algorithm/copy_if.hpp>

#include "MeshLib/Mesh.h"

namespace MeshGeoToolsLib
{
MeshNodesOnPoint::MeshNodesOnPoint(MeshLib::Mesh const& mesh,
                                   GeoLib::Grid<MeshLib::Node> const& mesh_grid,
                                   GeoLib::Point const& pnt,
                                   double epsilon_radius,
                                   SearchAllNodes search_all_nodes)
    : _mesh(mesh), _pnt(pnt)
{
    std::vector<std::size_t> vec_ids(
        mesh_grid.getPointsInEpsilonEnvironment(pnt, epsilon_radius));
    if (search_all_nodes == SearchAllNodes::Yes)
    {
        _msh_node_ids = vec_ids;
    }
    else
    {
        auto is_base_node = [this](std::size_t const id)
        {
            return MeshLib::isBaseNode(*_mesh.getNode(id),
                                       _mesh.getElementsConnectedToNode(id));
        };

        ranges::copy_if(
            vec_ids, std::back_inserter(_msh_node_ids), is_base_node);
    }
}

}  // end namespace MeshGeoToolsLib
