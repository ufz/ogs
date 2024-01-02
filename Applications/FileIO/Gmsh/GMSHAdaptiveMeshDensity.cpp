/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \date 2012-03-05
 * \author Thomas Fischer
 */

#include "GMSHAdaptiveMeshDensity.h"

#include <list>

#include "BaseLib/Logging.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"
#ifndef NDEBUG
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polyline.h"
#endif
#include "GeoLib/QuadTree.h"
#include "MathLib/MathTools.h"

namespace FileIO
{
namespace GMSH
{
GMSHAdaptiveMeshDensity::GMSHAdaptiveMeshDensity(double pnt_density,
                                                 double station_density,
                                                 std::size_t max_pnts_per_leaf)
    : _pnt_density(pnt_density),
      _station_density(station_density),
      _max_pnts_per_leaf(max_pnts_per_leaf),
      _quad_tree(nullptr)
{
}

GMSHAdaptiveMeshDensity::~GMSHAdaptiveMeshDensity()
{
    delete _quad_tree;
}

void GMSHAdaptiveMeshDensity::initialize(
    std::vector<GeoLib::Point const*> const& pnts)
{
    // *** QuadTree - determining bounding box
    DBUG(
        "GMSHAdaptiveMeshDensity::init(): computing axis aligned bounding box "
        "(2D) for quadtree.");

    GeoLib::Point min(*pnts[0]);
    GeoLib::Point max(*pnts[0]);
    std::size_t n_pnts(pnts.size());
    for (std::size_t k(1); k < n_pnts; k++)
    {
        for (std::size_t j(0); j < 2; j++)
        {
            if ((*(pnts[k]))[j] < min[j])
            {
                min[j] = (*(pnts[k]))[j];
            }
        }
        for (std::size_t j(0); j < 2; j++)
        {
            if ((*(pnts[k]))[j] > max[j])
            {
                max[j] = (*(pnts[k]))[j];
            }
        }
    }
    min[2] = 0.0;
    max[2] = 0.0;
    DBUG("GMSHAdaptiveMeshDensity::init(): \tok");

    // *** QuadTree - create object
    DBUG("GMSHAdaptiveMeshDensity::init(): Creating quadtree.");
    _quad_tree =
        new GeoLib::QuadTree<GeoLib::Point>(min, max, _max_pnts_per_leaf);
    DBUG("GMSHAdaptiveMeshDensity::init(): \tok.");

    // *** QuadTree - insert points
    addPoints(pnts);
}

void GMSHAdaptiveMeshDensity::addPoints(
    std::vector<GeoLib::Point const*> const& pnts)
{
    // *** QuadTree - insert points
    const std::size_t n_pnts(pnts.size());
    DBUG(
        "GMSHAdaptiveMeshDensity::addPoints(): Inserting {:d} points into "
        "quadtree.",
        n_pnts);
    for (std::size_t k(0); k < n_pnts; k++)
    {
        _quad_tree->addPoint(pnts[k]);
    }
    DBUG("GMSHAdaptiveMeshDensity::addPoints(): \tok.");
    _quad_tree->balance();
}

double GMSHAdaptiveMeshDensity::getMeshDensityAtPoint(
    GeoLib::Point const* const pnt) const
{
    GeoLib::Point ll;
    GeoLib::Point ur;
    _quad_tree->getLeaf(*pnt, ll, ur);
    return _pnt_density * (ur[0] - ll[0]);
}

double GMSHAdaptiveMeshDensity::getMeshDensityAtStation(
    GeoLib::Point const* const pnt) const
{
    GeoLib::Point ll;
    GeoLib::Point ur;
    _quad_tree->getLeaf(*pnt, ll, ur);
    return _station_density * (ur[0] - ll[0]);
}

void GMSHAdaptiveMeshDensity::getSteinerPoints(
    std::vector<GeoLib::Point*>& pnts, std::size_t additional_levels) const
{
    // get Steiner points
    std::size_t max_depth(0);
    _quad_tree->getMaxDepth(max_depth);

    std::list<GeoLib::QuadTree<GeoLib::Point>*> leaf_list;
    _quad_tree->getLeafs(leaf_list);

    for (std::list<GeoLib::QuadTree<GeoLib::Point>*>::const_iterator it(
             leaf_list.begin());
         it != leaf_list.end();
         ++it)
    {
        if ((*it)->getPoints().empty())
        {
            // compute point from square
            GeoLib::Point ll;
            GeoLib::Point ur;
            (*it)->getSquarePoints(ll, ur);
            if ((*it)->getDepth() + additional_levels > max_depth)
            {
                additional_levels = max_depth - (*it)->getDepth();
            }
            const std::size_t n_pnts_per_quad_dim = static_cast<std::size_t>(1)
                                                    << additional_levels;
            const double delta((ur[0] - ll[0]) / (2 * n_pnts_per_quad_dim));
            for (std::size_t i(0); i < n_pnts_per_quad_dim; i++)
            {
                for (std::size_t j(0); j < n_pnts_per_quad_dim; j++)
                {
                    pnts.push_back(new GeoLib::Point(
                        ll[0] + (2 * i + 1) * delta,
                        ll[1] + (2 * j + 1) * delta, 0.0, pnts.size()));
                }
            }
        }
    }
}

#ifndef NDEBUG
std::string GMSHAdaptiveMeshDensity::getQuadTreeGeometry(
    GeoLib::GEOObjects& geo_objs) const
{
    std::list<GeoLib::QuadTree<GeoLib::Point>*> leaf_list;
    _quad_tree->getLeafs(leaf_list);

    std::string quad_tree_geo("QuadTree");
    {
        std::vector<GeoLib::Point*> points{};
        for (auto const leaf : leaf_list)
        {
            // fetch corner points from leaf
            GeoLib::Point ll;
            GeoLib::Point ur;
            leaf->getSquarePoints(ll, ur);
            std::size_t const pnt_offset(points.size());
            points.push_back(new GeoLib::Point(ll, pnt_offset));
            points.push_back(
                new GeoLib::Point(ur[0], ll[1], 0.0, pnt_offset + 1));
            points.push_back(new GeoLib::Point(ur, pnt_offset + 2));
            points.push_back(
                new GeoLib::Point(ll[0], ur[1], 0.0, pnt_offset + 3));
        }
        geo_objs.addPointVec(std::move(points), quad_tree_geo,
                             GeoLib::PointVec::NameIdMap{});
    }
    auto& points = geo_objs.getPointVecObj(quad_tree_geo)->getVector();

    std::vector<GeoLib::Polyline*> polylines{};
    for (std::size_t l = 0; l < leaf_list.size(); ++l)
    {
        auto* polyline = new GeoLib::Polyline(points);
        for (std::size_t p = 0; p < 4; ++p)
        {
            polyline->addPoint(4 * l + p);
        }
        polyline->closePolyline();
        polylines.push_back(polyline);
    }
    geo_objs.addPolylineVec(std::move(polylines), quad_tree_geo,
                            GeoLib::PolylineVec::NameIdMap{});

    return quad_tree_geo;
}
#endif
}  // namespace GMSH

}  // end namespace FileIO
