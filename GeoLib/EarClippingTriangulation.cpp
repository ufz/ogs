/**
 * \file
 * \author Thomas Fischer
 * \date   2011-02-23
 * \brief  Implementation of the EarClippingTriangulation class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EarClippingTriangulation.h"

#include "BaseLib/uniqueInsert.h"

#include "MathLib/GeometricBasics.h"

#include "Polygon.h"
#include "Triangle.h"
#include "Point.h"

namespace GeoLib
{
EarClippingTriangulation::EarClippingTriangulation(const GeoLib::Polygon* polygon,
        std::list<GeoLib::Triangle> &triangles, bool rot)
{
    copyPolygonPoints (polygon);

    if (rot) {
        rotatePointsToXY (_pnts);
        ensureCWOrientation ();
    }

    initVertexList ();
    initLists ();
    clipEars ();

    std::vector<GeoLib::Point*> const& ref_pnts_vec (polygon->getPointsVec());
    std::list<GeoLib::Triangle>::const_iterator it (_triangles.begin());
    if (_original_orient == GeoLib::CW) {
        while (it != _triangles.end()) {
            const std::size_t i0 (polygon->getPointID ((*it)[0]));
            const std::size_t i1 (polygon->getPointID ((*it)[1]));
            const std::size_t i2 (polygon->getPointID ((*it)[2]));
            triangles.push_back (GeoLib::Triangle (ref_pnts_vec, i0, i1, i2));
            ++it;
        }
    } else {
        std::size_t n_pnts (_pnts.size() - 1);
        while (it != _triangles.end()) {
            const std::size_t i0 (polygon->getPointID (n_pnts - (*it)[0]));
            const std::size_t i1 (polygon->getPointID (n_pnts - (*it)[1]));
            const std::size_t i2 (polygon->getPointID (n_pnts - (*it)[2]));
            triangles.push_back (GeoLib::Triangle (ref_pnts_vec, i0, i1, i2));
            ++it;
        }
    }
}

EarClippingTriangulation::~EarClippingTriangulation()
{
    const std::size_t n_pnts (_pnts.size());
    for (std::size_t k(0); k<n_pnts; k++) {
        delete _pnts[k];
    }
}

void EarClippingTriangulation::copyPolygonPoints (const GeoLib::Polygon* polygon)
{
    // copy points - last point is identical to the first
    std::size_t n_pnts (polygon->getNumberOfPoints() - 1);
    for (std::size_t k(0); k < n_pnts; k++) {
        _pnts.push_back (new GeoLib::Point (*(polygon->getPoint(k))));
    }
}

void EarClippingTriangulation::ensureCWOrientation ()
{
    std::size_t n_pnts (_pnts.size());
    // get the left most upper point
    std::size_t min_x_max_y_idx (0); // for orientation check
    for (std::size_t k(0); k<n_pnts; k++) {
        if ((*(_pnts[k]))[0] <= (*(_pnts[min_x_max_y_idx]))[0]) {
            if ((*(_pnts[k]))[0] < (*(_pnts[min_x_max_y_idx]))[0]) {
                min_x_max_y_idx = k;
            } else {
                if ((*(_pnts[k]))[1] > (*(_pnts[min_x_max_y_idx]))[1]) {
                    min_x_max_y_idx = k;
                }
            }
        }
    }
    // determine orientation
    if (0 < min_x_max_y_idx && min_x_max_y_idx < n_pnts-1) {
        _original_orient = GeoLib::getOrientation (
            _pnts[min_x_max_y_idx-1], _pnts[min_x_max_y_idx], _pnts[min_x_max_y_idx+1]);
    } else {
        if (0 == min_x_max_y_idx) {
            _original_orient = GeoLib::getOrientation (_pnts[n_pnts-1], _pnts[0], _pnts[1]);
        } else {
            _original_orient = GeoLib::getOrientation (_pnts[n_pnts-2], _pnts[n_pnts-1], _pnts[0]);
        }
    }
    if (_original_orient == GeoLib::CCW) {
        // switch orientation
        for (std::size_t k(0); k<n_pnts/2; k++) {
            std::swap (_pnts[k], _pnts[n_pnts-1-k]);
        }
    }
}

bool EarClippingTriangulation::isEar(std::size_t v0, std::size_t v1, std::size_t v2) const
{
    for (unsigned long it : _vertex_list)
    {
        if (it != v0 && it != v1 && it != v2)
        {
            if (MathLib::isPointInTriangle(*_pnts[it], *_pnts[v0], *_pnts[v1],
                                           *_pnts[v2]))
            {
                return false;
            }
        }
    }
    return true;
}

void EarClippingTriangulation::initVertexList ()
{
    std::size_t n_pnts (_pnts.size());
    for (std::size_t k(0); k<n_pnts; k++)
        _vertex_list.push_back (k);
}

void EarClippingTriangulation::initLists ()
{
    // go through points checking ccw, cw or collinear order and identifying ears
    std::list<std::size_t>::iterator it (_vertex_list.begin()), prev(_vertex_list.end()), next;
    --prev;
    next = it;
    ++next;
    GeoLib::Orientation orientation;
    bool first_run(true); // saves special handling of the last case with identical code
    while (_vertex_list.size() >= 3 && first_run) {
        if (next == _vertex_list.end()) {
            first_run = false;
            next = _vertex_list.begin();
        }
        orientation  = getOrientation (_pnts[*prev], _pnts[*it], _pnts[*next]);
        if (orientation == GeoLib::COLLINEAR) {
            WARN("EarClippingTriangulation::initLists(): collinear points (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)",
                    (*_pnts[*prev])[0], (*_pnts[*prev])[1], (*_pnts[*prev])[2],
                    (*_pnts[*it])[0], (*_pnts[*it])[1], (*_pnts[*it])[2],
                    (*_pnts[*next])[0], (*_pnts[*next])[1], (*_pnts[*next])[2]);
            it = _vertex_list.erase (it);
            ++next;
        } else {
            if (orientation == GeoLib::CW) {
                _convex_vertex_list.push_back (*it);
                if (isEar (*prev, *it, *next))
                    _ear_list.push_back (*it);
            }
            prev = it;
            it = next;
            ++next;
        }
    }
}

void EarClippingTriangulation::clipEars()
{
    std::list<std::size_t>::iterator it, prev, next;
    // *** clip an ear
    while (_vertex_list.size() > 3) {
        // pop ear from list
        std::size_t ear = _ear_list.front();
        _ear_list.pop_front();
        // remove ear tip from _convex_vertex_list
        _convex_vertex_list.remove(ear);

        // remove ear from vertex_list, apply changes to _ear_list, _convex_vertex_list
        bool nfound(true);
        it = _vertex_list.begin();
        prev = _vertex_list.end();
        --prev;
        while (nfound && it != _vertex_list.end()) {
            if (*it == ear) {
                nfound = false;
                it = _vertex_list.erase(it); // remove ear tip
                next = it;
                if (next == _vertex_list.end()) {
                    next = _vertex_list.begin();
                    prev = _vertex_list.end();
                    --prev;
                }
                // add triangle
                _triangles.push_back(GeoLib::Triangle(_pnts, *prev, *next, ear));

                // check the orientation of prevprev, prev, next
                std::list<std::size_t>::iterator prevprev;
                if (prev == _vertex_list.begin()) {
                    prevprev = _vertex_list.end();
                } else {
                    prevprev = prev;
                }
                --prevprev;

                // apply changes to _convex_vertex_list and _ear_list looking "backward"
                GeoLib::Orientation orientation = GeoLib::getOrientation(_pnts[*prevprev], _pnts[*prev],
                        _pnts[*next]);
                if (orientation == GeoLib::CW) {
                    BaseLib::uniquePushBack(_convex_vertex_list, *prev);
                    // prev is convex
                    if (isEar(*prevprev, *prev, *next)) {
                        // prev is an ear tip
                        BaseLib::uniquePushBack(_ear_list, *prev);
                    } else {
                        // if necessary remove prev
                        _ear_list.remove(*prev);
                    }
                } else {
                    // prev is not convex => reflex or collinear
                    _convex_vertex_list.remove(*prev);
                    _ear_list.remove(*prev);
                    if (orientation == GeoLib::COLLINEAR) {
                        prev = _vertex_list.erase(prev);
                        if (prev == _vertex_list.begin()) {
                            prev = _vertex_list.end();
                            --prev;
                        } else {
                            --prev;
                        }
                    }
                }

                // check the orientation of prev, next, nextnext
                std::list<std::size_t>::iterator nextnext,
                        help_it(_vertex_list.end());
                --help_it;
                if (next == help_it) {
                    nextnext = _vertex_list.begin();
                } else {
                    nextnext = next;
                    ++nextnext;
                }

                // apply changes to _convex_vertex_list and _ear_list looking "forward"
                orientation = getOrientation(_pnts[*prev], _pnts[*next],
                        _pnts[*nextnext]);
                if (orientation == GeoLib::CW) {
                    BaseLib::uniquePushBack(_convex_vertex_list, *next);
                    // next is convex
                    if (isEar(*prev, *next, *nextnext)) {
                        // next is an ear tip
                        BaseLib::uniquePushBack(_ear_list, *next);
                    } else {
                        // if necessary remove *next
                        _ear_list.remove(*next);
                    }
                } else {
                    // next is not convex => reflex or collinear
                    _convex_vertex_list.remove(*next);
                    _ear_list.remove(*next);
                    if (orientation == GeoLib::COLLINEAR) {
                        next = _vertex_list.erase(next);
                        if (next == _vertex_list.end())
                            next = _vertex_list.begin();
                    }
                }
            } else {
                prev = it;
                ++it;
            }
        }

    }

    // add last triangle
    next = _vertex_list.begin();
    prev = next;
    ++next;
    if (next == _vertex_list.end())
        return;
    it = next;
    ++next;
    if (next == _vertex_list.end())
        return;

    if (getOrientation(_pnts[*prev], _pnts[*it], _pnts[*next]) == GeoLib::CCW)
        _triangles.push_back(GeoLib::Triangle(_pnts, *prev, *it, *next));
    else
        _triangles.push_back(GeoLib::Triangle(_pnts, *prev, *next, *it));
}

} // end namespace GeoLib
