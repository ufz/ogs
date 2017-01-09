/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "MeshNodesOnPoint.h"

#include "MeshLib/Mesh.h"

namespace MeshGeoToolsLib
{

MeshNodesOnPoint::MeshNodesOnPoint(MeshLib::Mesh const& mesh, GeoLib::Grid<MeshLib::Node> const &mesh_grid,
        GeoLib::Point const& pnt, double epsilon_radius, bool search_all_nodes)
: _mesh(mesh), _pnt(pnt)
{
    std::vector<std::size_t> vec_ids(mesh_grid.getPointsInEpsilonEnvironment(pnt, epsilon_radius));
    if (search_all_nodes)
        _msh_node_ids = vec_ids;
    else {
        for (auto id : vec_ids)
            if (mesh.isBaseNode(id))
                _msh_node_ids.push_back(id);
    }
}

} // end namespace MeshGeoToolsLib

