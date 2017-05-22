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

#include <typeinfo>

#include "HeuristicSearchLength.h"
#include "MeshNodesAlongPolyline.h"
#include "MeshNodesAlongSurface.h"
#include "MeshNodesOnPoint.h"

#include <logog/include/logog.hpp>

#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{

std::vector<std::unique_ptr<MeshNodeSearcher>> MeshNodeSearcher::_mesh_node_searchers;

MeshNodeSearcher::MeshNodeSearcher(
    MeshLib::Mesh const& mesh,
    std::unique_ptr<MeshGeoToolsLib::SearchLength>&& search_length_algorithm,
    SearchAllNodes search_all_nodes)
    : _mesh(mesh),
      _mesh_grid(_mesh.getNodes().cbegin(), _mesh.getNodes().cend()),
      _search_length_algorithm(std::move(search_length_algorithm)),
      _search_all_nodes(search_all_nodes)
{
    DBUG("The search length for mesh \"%s\" is %e.",
        _mesh.getName().c_str(), _search_length_algorithm->getSearchLength());
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

std::vector<std::size_t> MeshNodeSearcher::getMeshNodeIDs(
    GeoLib::GeoObject const& geoObj) const
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

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsForPoint(
    GeoLib::Point const& pnt) const
{
    return getMeshNodesOnPoint(pnt).getNodeIDs();
}

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsAlongPolyline(
    GeoLib::Polyline const& ply) const
{
    return getMeshNodesAlongPolyline(ply).getNodeIDs();
}

std::vector<std::size_t> const& MeshNodeSearcher::getMeshNodeIDsAlongSurface(
    GeoLib::Surface const& sfc) const
{
    return getMeshNodesAlongSurface(sfc).getNodeIDs();
}

MeshNodesOnPoint& MeshNodeSearcher::getMeshNodesOnPoint(
    GeoLib::Point const& pnt) const
{
    std::vector<MeshNodesOnPoint*>::const_iterator it(_mesh_nodes_on_points.begin());
    for (; it != _mesh_nodes_on_points.end(); ++it) {
        if (&(*it)->getPoint() == &pnt) {
            return *(*it);
        }
    }

    _mesh_nodes_on_points.push_back(
        new MeshNodesOnPoint(_mesh,
                             _mesh_grid,
                             pnt,
                             _search_length_algorithm->getSearchLength(),
                             _search_all_nodes));
    return *_mesh_nodes_on_points.back();
}

MeshNodesAlongPolyline& MeshNodeSearcher::getMeshNodesAlongPolyline(
    GeoLib::Polyline const& ply) const
{
    std::vector<MeshNodesAlongPolyline*>::const_iterator it(_mesh_nodes_along_polylines.begin());
    for (; it != _mesh_nodes_along_polylines.end(); ++it) {
        if (&(*it)->getPolyline() == &ply) {
            // we calculated mesh nodes for this polyline already
            return *(*it);
        }
    }

    // compute nodes (and supporting points) along polyline
    _mesh_nodes_along_polylines.push_back(new MeshNodesAlongPolyline(
        _mesh, ply, _search_length_algorithm->getSearchLength(),
        _search_all_nodes));
    return *_mesh_nodes_along_polylines.back();
}

MeshNodesAlongSurface& MeshNodeSearcher::getMeshNodesAlongSurface(
    GeoLib::Surface const& sfc) const
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
        new MeshNodesAlongSurface(_mesh,
                                  sfc,
                                  _search_length_algorithm->getSearchLength(),
                                  _search_all_nodes));
    return *_mesh_nodes_along_surfaces.back();
}

MeshNodeSearcher const& MeshNodeSearcher::getMeshNodeSearcher(
    MeshLib::Mesh const& mesh,
    std::unique_ptr<MeshGeoToolsLib::SearchLength>&& search_length_algorithm)
{
    std::size_t const mesh_id = mesh.getID();
    if (_mesh_node_searchers.size() < mesh_id+1)
        _mesh_node_searchers.resize(mesh_id+1);

    if (_mesh_node_searchers[mesh_id])
    {
        auto const& m = *_mesh_node_searchers[mesh_id];
        // return searcher if search length algorithm and the returned search
        // lenght are the same, else recreate the searcher
        if (typeid(m._search_length_algorithm) ==
                typeid(search_length_algorithm) &&
            m._search_length_algorithm->getSearchLength() ==
                search_length_algorithm->getSearchLength())
        {
            return m;
        }
    }

    _mesh_node_searchers[mesh_id] =
        std::make_unique<MeshGeoToolsLib::MeshNodeSearcher>(
            mesh, std::move(search_length_algorithm), SearchAllNodes::Yes);

    return *_mesh_node_searchers[mesh_id];
}

std::size_t
MeshNodeSearcher::getMeshId() const
{
    return _mesh.getID();
}

} // end namespace MeshGeoToolsLib
