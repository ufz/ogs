/**
 * \file
 * \date Oct 24, 2013
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshNodeSearcher.h"

#include <typeinfo>
#include <sstream>

#include "HeuristicSearchLength.h"
#include "MeshNodesAlongPolyline.h"
#include "MeshNodesAlongSurface.h"
#include "MeshNodesOnPoint.h"

#include "BaseLib/Logging.h"

#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
std::vector<std::unique_ptr<MeshNodeSearcher>>
    MeshNodeSearcher::_mesh_node_searchers;

MeshNodeSearcher::MeshNodeSearcher(
    MeshLib::Mesh const& mesh,
    std::unique_ptr<MeshGeoToolsLib::SearchLength>&& search_length_algorithm,
    SearchAllNodes search_all_nodes)
    : _mesh(mesh),
      _mesh_grid(_mesh.getNodes().cbegin(), _mesh.getNodes().cend()),
      _search_length_algorithm(std::move(search_length_algorithm)),
      _search_all_nodes(search_all_nodes)
{
    DBUG("The search length for mesh '{:s}' is {:e}.", _mesh.getName(),
         _search_length_algorithm->getSearchLength());
}

MeshNodeSearcher::~MeshNodeSearcher()
{
    for (auto pointer : _mesh_nodes_on_points)
    {
        delete pointer;
    }
    for (auto pointer : _mesh_nodes_along_polylines)
    {
        delete pointer;
    }
    for (auto pointer : _mesh_nodes_along_surfaces)
    {
        delete pointer;
    }
}

template <typename CacheType, typename GeometryType>
std::vector<std::size_t> const& getMeshNodeIDs(
    std::vector<CacheType*>& cached_elements,
    std::function<GeometryType(CacheType const&)> getCachedItem,
    GeometryType const& item, MeshLib::Mesh const& mesh,
    GeoLib::Grid<MeshLib::Node> const& mesh_grid,
    double const search_length,
    SearchAllNodes const search_all_nodes)
{
    if (auto const it = find_if(cbegin(cached_elements), cend(cached_elements),
                                [&](auto const& element) {
                                    return getCachedItem(*element) == item;
                                });
        it != cend(cached_elements))
    {
        return (*it)->getNodeIDs();
    }
    // search IDs for geometry object
    if constexpr (std::is_convertible<GeometryType, GeoLib::Point>::value)
    {
        cached_elements.push_back(new CacheType(
            mesh, mesh_grid, item, search_length, search_all_nodes));
    }
    else
    {
        cached_elements.push_back(
            new CacheType(mesh, item, search_length, search_all_nodes));
    }
    return cached_elements.back()->getNodeIDs();
}

std::vector<std::size_t> MeshNodeSearcher::getMeshNodeIDs(
    GeoLib::GeoObject const& geoObj) const
{
    std::vector<std::size_t> vec_nodes;
    switch (geoObj.getGeoType())
    {
        case GeoLib::GEOTYPE::POINT:
        {
            std::function<GeoLib::Point(MeshNodesOnPoint const&)>
                get_cached_item_function = &MeshNodesOnPoint::getPoint;
            return MeshGeoToolsLib::getMeshNodeIDs(
                _mesh_nodes_on_points, get_cached_item_function,
                *static_cast<const GeoLib::Point*>(&geoObj), _mesh, _mesh_grid,
                _search_length_algorithm->getSearchLength(), _search_all_nodes);
        }
        case GeoLib::GEOTYPE::POLYLINE:
        {
            std::function<GeoLib::Polyline(MeshNodesAlongPolyline const&)>
                get_cached_item_function = &MeshNodesAlongPolyline::getPolyline;
            return MeshGeoToolsLib::getMeshNodeIDs(
                _mesh_nodes_along_polylines, get_cached_item_function,
                *static_cast<const GeoLib::Polyline*>(&geoObj), _mesh,
                _mesh_grid, _search_length_algorithm->getSearchLength(),
                _search_all_nodes);
        }
        case GeoLib::GEOTYPE::SURFACE:
        {
            std::function<GeoLib::Surface(MeshNodesAlongSurface const&)>
                get_cached_item_function = &MeshNodesAlongSurface::getSurface;
            return MeshGeoToolsLib::getMeshNodeIDs(
                _mesh_nodes_along_surfaces, get_cached_item_function,
                *static_cast<const GeoLib::Surface*>(&geoObj), _mesh,
                _mesh_grid, _search_length_algorithm->getSearchLength(),
                _search_all_nodes);
        }
        default:
            break;
    }
    return vec_nodes;
}

std::vector<std::size_t> MeshNodeSearcher::getMeshNodeIDs(
    std::vector<MathLib::Point3dWithID*> const& points) const
{
    double const epsilon_radius = _search_length_algorithm->getSearchLength();

    std::vector<std::size_t> node_ids;
    node_ids.reserve(points.size());

    for (auto const* const p_ptr : points)
    {
        auto const& p = *p_ptr;
        std::vector<std::size_t> const ids =
            _mesh_grid.getPointsInEpsilonEnvironment(p, epsilon_radius);
        if (ids.empty())
        {
            OGS_FATAL(
                "No nodes could be found in the mesh for point {:d} : ({:g}, "
                "{:g}, {:g}) in {:g} epsilon radius in the mesh '{:s}'",
                p.getID(), p[0], p[1], p[2], epsilon_radius, _mesh.getName());
        }
        if (ids.size() != 1)
        {
            std::stringstream ss;
            auto const& bulk_nodes = _mesh.getNodes();
            for (auto const id : ids)
            {
                ss << "- bulk node: " << (*bulk_nodes[id])[0] << ", "
                   << (*bulk_nodes[id])[1] << ", " << (*bulk_nodes[id])[2]
                   << ", distance: "
                   << std::sqrt(MathLib::sqrDist(*bulk_nodes[id], p)) << "\n";
            }
            OGS_FATAL(
                "Found {:d} nodes in the mesh for point {:d} : ({:g}, {:g}, "
                "{:g}) in {:g} epsilon radius in the mesh '{:s}'. Expected to "
                "find exactly one node.\n{:s}",
                ids.size(), p.getID(), p[0], p[1], p[2], epsilon_radius,
                _mesh.getName(), ss.str());
        }
        node_ids.push_back(ids.front());
    }
    return node_ids;
}

MeshNodeSearcher const& MeshNodeSearcher::getMeshNodeSearcher(
    MeshLib::Mesh const& mesh,
    std::unique_ptr<MeshGeoToolsLib::SearchLength>&& search_length_algorithm)
{
    std::size_t const mesh_id = mesh.getID();
    if (_mesh_node_searchers.size() < mesh_id + 1)
    {
        _mesh_node_searchers.resize(mesh_id + 1);
    }

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

}  // end namespace MeshGeoToolsLib
