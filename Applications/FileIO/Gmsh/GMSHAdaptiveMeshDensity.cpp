/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file GMSHAdaptiveMeshDensity.cpp
 * @date 2012-03-05
 * @author Thomas Fischer
 */

#include "GMSHAdaptiveMeshDensity.h"

#include <list>

#include <logog/include/logog.hpp>

#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"
#ifndef NDEBUG
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

void GMSHAdaptiveMeshDensity::initialize(std::vector<GeoLib::Point const*> const& pnts)
{
    // *** QuadTree - determining bounding box
    DBUG("GMSHAdaptiveMeshDensity::init(): computing axis aligned bounding box (2D) for quadtree.");

    GeoLib::Point min(*pnts[0]), max(*pnts[0]);
    std::size_t n_pnts(pnts.size());
    for (std::size_t k(1); k<n_pnts; k++) {
        for (std::size_t j(0); j < 2; j++)
            if ((*(pnts[k]))[j] < min[j])
                min[j] = (*(pnts[k]))[j];
        for (std::size_t j(0); j < 2; j++)
            if ((*(pnts[k]))[j] > max[j])
                max[j] = (*(pnts[k]))[j];
    }
    min[2] = 0.0;
    max[2] = 0.0;
    DBUG("GMSHAdaptiveMeshDensity::init(): \tok");

    // *** QuadTree - create object
    DBUG("GMSHAdaptiveMeshDensity::init(): Creating quadtree.");
    _quad_tree = new GeoLib::QuadTree<GeoLib::Point> (min, max, _max_pnts_per_leaf);
    DBUG("GMSHAdaptiveMeshDensity::init(): \tok.");

    // *** QuadTree - insert points
    addPoints(pnts);
}

void GMSHAdaptiveMeshDensity::addPoints(std::vector<GeoLib::Point const*> const& pnts)
{
    // *** QuadTree - insert points
    const std::size_t n_pnts(pnts.size());
    DBUG("GMSHAdaptiveMeshDensity::addPoints(): Inserting %d points into quadtree.", n_pnts);
    for (std::size_t k(0); k < n_pnts; k++)
        _quad_tree->addPoint(pnts[k]);
    DBUG("GMSHAdaptiveMeshDensity::addPoints(): \tok.");
    _quad_tree->balance();
}

double GMSHAdaptiveMeshDensity::getMeshDensityAtPoint(GeoLib::Point const* const pnt) const
{
    GeoLib::Point ll, ur;
    _quad_tree->getLeaf(*pnt, ll, ur);
    return _pnt_density * (ur[0] - ll[0]);
}

double GMSHAdaptiveMeshDensity::getMeshDensityAtStation(GeoLib::Point const* const pnt) const
{
    GeoLib::Point ll, ur;
    _quad_tree->getLeaf(*pnt, ll, ur);
    return _station_density * (ur[0] - ll[0]);
}

void GMSHAdaptiveMeshDensity::getSteinerPoints (std::vector<GeoLib::Point*> & pnts,
                                                std::size_t additional_levels) const
{
    // get Steiner points
    std::size_t max_depth(0);
    _quad_tree->getMaxDepth(max_depth);

    std::list<GeoLib::QuadTree<GeoLib::Point>*> leaf_list;
    _quad_tree->getLeafs(leaf_list);

    for (std::list<GeoLib::QuadTree<GeoLib::Point>*>::const_iterator it(leaf_list.begin()); it
                    != leaf_list.end(); ++it) {
        if ((*it)->getPoints().empty()) {
            // compute point from square
            GeoLib::Point ll, ur;
            (*it)->getSquarePoints(ll, ur);
            if ((*it)->getDepth() + additional_levels > max_depth) {
                additional_levels = max_depth - (*it)->getDepth();
            }
            const std::size_t n_pnts_per_quad_dim = static_cast<std::size_t>(1)
                                                    << additional_levels;
            const double delta ((ur[0] - ll[0]) / (2 * n_pnts_per_quad_dim));
            for (std::size_t i(0); i<n_pnts_per_quad_dim; i++) {
                for (std::size_t j(0); j<n_pnts_per_quad_dim; j++) {
                    pnts.push_back(new GeoLib::Point(
                        ll[0] + (2 * i + 1) * delta,
                        ll[1] + (2 * j + 1) * delta, 0.0, pnts.size()));
                }
            }

        }
    }
}

#ifndef NDEBUG
void GMSHAdaptiveMeshDensity::getQuadTreeGeometry(std::vector<GeoLib::Point*> &pnts,
                                                  std::vector<GeoLib::Polyline*> &plys) const
{
    std::list<GeoLib::QuadTree<GeoLib::Point>*> leaf_list;
    _quad_tree->getLeafs(leaf_list);

    for (std::list<GeoLib::QuadTree<GeoLib::Point>*>::const_iterator it(leaf_list.begin()); it
        != leaf_list.end(); ++it) {
        // fetch corner points from leaf
        GeoLib::Point ll, ur;
        (*it)->getSquarePoints(ll, ur);
        std::size_t const pnt_offset (pnts.size());
        pnts.push_back(new GeoLib::Point(ll, pnt_offset));
        pnts.push_back(new GeoLib::Point(ur[0], ll[1], 0.0, pnt_offset+1));
        pnts.push_back(new GeoLib::Point(ur, pnt_offset+2));
        pnts.push_back(new GeoLib::Point(ll[0], ur[1], 0.0, pnt_offset+3));
        plys.push_back(new GeoLib::Polyline(pnts));
        plys[plys.size() - 1]->addPoint(pnt_offset);
        plys[plys.size() - 1]->addPoint(pnt_offset + 1);
        plys[plys.size() - 1]->addPoint(pnt_offset + 2);
        plys[plys.size() - 1]->addPoint(pnt_offset + 3);
        plys[plys.size() - 1]->addPoint(pnt_offset);
    }
}
#endif
}

} // end namespace FileIO
