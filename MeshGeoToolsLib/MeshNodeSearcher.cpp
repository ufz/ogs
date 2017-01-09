/**
 * @date Oct 24, 2013
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "MeshNodeSearcher.h"
#include "HeuristicSearchLength.h"
#include "MeshNodesAlongPolyline.h"
#include "MeshNodesAlongSurface.h"
#include "MeshNodesOnPoint.h"

#include <logog/include/logog.hpp>

// GeoLib
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"

// MeshLib
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{

std::vector<std::unique_ptr<MeshNodeSearcher>> MeshNodeSearcher::_mesh_node_searchers;


MeshNodeSearcher::MeshNodeSearcher(MeshLib::Mesh const& mesh,
    MeshGeoToolsLib::SearchLength const& search_length_algorithm, bool search_all_nodes) :
        _mesh(mesh), _mesh_grid(_mesh.getNodes().cbegin(), _mesh.getNodes().cend()),
        _search_length(0.0), _search_all_nodes(search_all_nodes)
{
    DBUG("Constructing MeshNodeSearcher obj.");
    _search_length = search_length_algorithm.getSearchLength();

    INFO("Calculated search length for mesh \"%s\" is %e.",
        _mesh.getName().c_str(), _search_length);
}

MeshNodeSearcher::~MeshNodeSearcher()
{
    for (auto pointer : _mesh_nodes_on_points)
        delete pointer;
    for (auto pointer : _mesh_nodes_along_polylines)
        delete pointer;
    for (auto pointer : _mesh_nodes_along_surfaces)
        delete pointer;
}

std::vector<std::size_t> MeshNodeSearcher::getMeshNodeIDs(GeoLib::GeoObject const& geoObj)
{
    std::vector<std::size_t> vec_nodes;
    switch (geoObj.getGeoType()) {
    case GeoLib::GEOTYPE::POINT:
    {
        vec_nodes = this->getMeshNodeIDsForPoint(*static_cast<const GeoLib::Point*>(&geoObj));
        break;
    }
    case GeoLib::GEOTYPE::POLYLINE:
        vec_nodes = this->getMeshNodeIDsAlongPolyline(*static_cast<const GeoLib::Polyline*>(&geoObj));
        break;
    case GeoLib::GEOTYPE::SURFACE:
        vec_nodes = this->getMeshNodeIDsAlongSurface(*static_cast<const GeoLib::Surface*>(&geoObj));
        break;
    default:
        break;
    }
    return vec_nodes;
}

std::vector<std::size_t> const&
MeshNodeSearcher::getMeshNodeIDsForPoint(GeoLib::Point const& pnt)
{
    return getMeshNodesOnPoint(pnt).getNodeIDs();
}

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsAlongPolyline(
        GeoLib::Polyline const& ply)
{
    return getMeshNodesAlongPolyline(ply).getNodeIDs();
}

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsAlongSurface(GeoLib::Surface const& sfc)
{
    return getMeshNodesAlongSurface(sfc).getNodeIDs();
}

MeshNodesOnPoint& MeshNodeSearcher::getMeshNodesOnPoint(GeoLib::Point const& pnt)
{
    std::vector<MeshNodesOnPoint*>::const_iterator it(_mesh_nodes_on_points.begin());
    for (; it != _mesh_nodes_on_points.end(); ++it) {
        if (&(*it)->getPoint() == &pnt) {
            return *(*it);
        }
    }

    _mesh_nodes_on_points.push_back(
            new MeshNodesOnPoint(_mesh, _mesh_grid, pnt, _search_length, _search_all_nodes));
    return *_mesh_nodes_on_points.back();
}

MeshNodesAlongPolyline& MeshNodeSearcher::getMeshNodesAlongPolyline(GeoLib::Polyline const& ply)
{
    std::vector<MeshNodesAlongPolyline*>::const_iterator it(_mesh_nodes_along_polylines.begin());
    for (; it != _mesh_nodes_along_polylines.end(); ++it) {
        if (&(*it)->getPolyline() == &ply) {
            // we calculated mesh nodes for this polyline already
            return *(*it);
        }
    }

    // compute nodes (and supporting points) along polyline
    _mesh_nodes_along_polylines.push_back(
            new MeshNodesAlongPolyline(_mesh, ply, _search_length, _search_all_nodes));
    return *_mesh_nodes_along_polylines.back();
}

MeshNodesAlongSurface& MeshNodeSearcher::getMeshNodesAlongSurface(GeoLib::Surface const& sfc)
{
    std::vector<MeshNodesAlongSurface*>::const_iterator it(_mesh_nodes_along_surfaces.begin());
    for (; it != _mesh_nodes_along_surfaces.end(); ++it) {
        if (&(*it)->getSurface() == &sfc) {
            // we calculated mesh nodes for this polyline already
            return *(*it);
        }
    }

    // compute nodes (and supporting points) along polyline
    _mesh_nodes_along_surfaces.push_back(
            new MeshNodesAlongSurface(_mesh, sfc, _search_length, _search_all_nodes));
    return *_mesh_nodes_along_surfaces.back();
}

MeshNodeSearcher&
MeshNodeSearcher::getMeshNodeSearcher(MeshLib::Mesh const& mesh)
{
    std::size_t const mesh_id = mesh.getID();
    if (_mesh_node_searchers.size() < mesh_id+1)
        _mesh_node_searchers.resize(mesh_id+1);

    if (!_mesh_node_searchers[mesh_id])
        _mesh_node_searchers[mesh_id].reset(
            new MeshGeoToolsLib::MeshNodeSearcher(mesh));

    return *_mesh_node_searchers[mesh_id];
}

std::size_t
MeshNodeSearcher::getMeshId() const
{
    return _mesh.getID();
}

} // end namespace MeshGeoToolsLib
