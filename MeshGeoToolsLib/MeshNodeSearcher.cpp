/**
 * \file
 * \date Oct 24, 2013
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    MeshNodeSearcher::mesh_node_searchers_;

MeshNodeSearcher::MeshNodeSearcher(
    MeshLib::Mesh const& mesh,
    std::unique_ptr<MeshGeoToolsLib::SearchLength>&& search_length_algorithm,
    SearchAllNodes search_all_nodes)
    : mesh_(mesh),
      mesh_grid_(mesh_.getNodes().cbegin(), mesh_.getNodes().cend()),
      search_length_algorithm_(std::move(search_length_algorithm)),
      search_all_nodes_(search_all_nodes)
{
    DBUG("The search length for mesh '{:s}' is {:e}.", mesh_.getName(),
         search_length_algorithm_->getSearchLength());
}

MeshNodeSearcher::~MeshNodeSearcher()
{
    for (auto pointer : mesh_nodes_on_points_)
    {
        delete pointer;
    }
    for (auto pointer : mesh_nodes_along_polylines_)
    {
        delete pointer;
    }
    for (auto pointer : mesh_nodes_along_surfaces_)
    {
        delete pointer;
    }
}

std::vector<std::size_t> MeshNodeSearcher::getMeshNodeIDs(
    GeoLib::GeoObject const& geoObj) const
{
    std::vector<std::size_t> vec_nodes;
    switch (geoObj.getGeoType())
    {
        case GeoLib::GEOTYPE::POINT:
        {
            vec_nodes = this->getMeshNodeIDsForPoint(
                *static_cast<const GeoLib::Point*>(&geoObj));
            break;
        }
        case GeoLib::GEOTYPE::POLYLINE:
            vec_nodes = this->getMeshNodeIDsAlongPolyline(
                *static_cast<const GeoLib::Polyline*>(&geoObj));
            break;
        case GeoLib::GEOTYPE::SURFACE:
            vec_nodes = this->getMeshNodeIDsAlongSurface(
                *static_cast<const GeoLib::Surface*>(&geoObj));
            break;
        default:
            break;
    }
    return vec_nodes;
}

std::vector<std::size_t> MeshNodeSearcher::getMeshNodeIDs(
    std::vector<MathLib::Point3dWithID*> const& points) const
{
    double const epsilon_radius = search_length_algorithm_->getSearchLength();

    std::vector<std::size_t> node_ids;
    node_ids.reserve(points.size());

    for (auto const* const p_ptr : points)
    {
        auto const& p = *p_ptr;
        std::vector<std::size_t> const ids =
            mesh_grid_.getPointsInEpsilonEnvironment(p, epsilon_radius);
        if (ids.empty())
        {
            OGS_FATAL(
                "No nodes could be found in the mesh for point {:d} : ({:g}, "
                "{:g}, {:g}) in {:g} epsilon radius in the mesh '{:s}'",
                p.getID(), p[0], p[1], p[2], epsilon_radius, mesh_.getName());
        }
        if (ids.size() != 1)
        {
            std::stringstream ss;
            auto const& bulk_nodes = mesh_.getNodes();
            for (auto const id : ids)
            {
                ss << "- bulk node: " << (*bulk_nodes[id]) << ", distance: "
                   << std::sqrt(MathLib::sqrDist(bulk_nodes[id]->getCoords(),
                                                 p.getCoords()))
                   << "\n";
            }
            OGS_FATAL(
                "Found {:d} nodes in the mesh for point {:d} : ({:g}, {:g}, "
                "{:g}) in {:g} epsilon radius in the mesh '{:s}'. Expected to "
                "find exactly one node.\n{:s}",
                ids.size(), p.getID(), p[0], p[1], p[2], epsilon_radius,
                mesh_.getName(), ss.str());
        }
        node_ids.push_back(ids.front());
    }
    return node_ids;
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
    for (auto const& mesh_nodes_on_point : mesh_nodes_on_points_)
    {
        if (&(mesh_nodes_on_point->getPoint()) == &pnt)
        {
            return *mesh_nodes_on_point;
        }
    }

    mesh_nodes_on_points_.push_back(
        new MeshNodesOnPoint(mesh_,
                             mesh_grid_,
                             pnt,
                             search_length_algorithm_->getSearchLength(),
                             search_all_nodes_));
    return *mesh_nodes_on_points_.back();
}

MeshNodesAlongPolyline& MeshNodeSearcher::getMeshNodesAlongPolyline(
    GeoLib::Polyline const& ply) const
{
    for (auto const& mesh_nodes_along_polyline : mesh_nodes_along_polylines_)
    {
        if (&(mesh_nodes_along_polyline->getPolyline()) == &ply)
        {
            return *mesh_nodes_along_polyline;
        }
    }

    // compute nodes (and supporting points) along polyline
    mesh_nodes_along_polylines_.push_back(new MeshNodesAlongPolyline(
        mesh_, ply, search_length_algorithm_->getSearchLength(),
        search_all_nodes_));
    return *mesh_nodes_along_polylines_.back();
}

MeshNodesAlongSurface& MeshNodeSearcher::getMeshNodesAlongSurface(
    GeoLib::Surface const& sfc) const
{
    for (auto const& mesh_nodes_along_surface : mesh_nodes_along_surfaces_)
    {
        if (&(mesh_nodes_along_surface->getSurface()) == &sfc)
        {
            return *mesh_nodes_along_surface;
        }
    }

    // compute nodes (and supporting points) on surface
    mesh_nodes_along_surfaces_.push_back(
        new MeshNodesAlongSurface(mesh_,
                                  sfc,
                                  search_length_algorithm_->getSearchLength(),
                                  search_all_nodes_));
    return *mesh_nodes_along_surfaces_.back();
}

MeshNodeSearcher const& MeshNodeSearcher::getMeshNodeSearcher(
    MeshLib::Mesh const& mesh,
    std::unique_ptr<MeshGeoToolsLib::SearchLength>&& search_length_algorithm)
{
    std::size_t const mesh_id = mesh.getID();
    if (mesh_node_searchers_.size() < mesh_id + 1)
    {
        mesh_node_searchers_.resize(mesh_id + 1);
    }

    if (mesh_node_searchers_[mesh_id])
    {
        auto const& m = *mesh_node_searchers_[mesh_id];
        // return searcher if search length algorithm and the returned search
        // lenght are the same, else recreate the searcher
        if (typeid(m.search_length_algorithm_) ==
                typeid(search_length_algorithm) &&
            m.search_length_algorithm_->getSearchLength() ==
                search_length_algorithm->getSearchLength())
        {
            return m;
        }
    }

    mesh_node_searchers_[mesh_id] =
        std::make_unique<MeshGeoToolsLib::MeshNodeSearcher>(
            mesh, std::move(search_length_algorithm), SearchAllNodes::Yes);

    return *mesh_node_searchers_[mesh_id];
}

}  // end namespace MeshGeoToolsLib
