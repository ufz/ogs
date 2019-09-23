/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <list>

#include <logog/include/logog.hpp>

#include "Surface.h"

// GeoLib
#include "AABB.h"
#include "Polygon.h"
#include "SurfaceGrid.h"
#include "AnalyticalGeometry.h"

#include "Triangle.h"
#include "Polyline.h"

namespace GeoLib
{
Surface::Surface(const std::vector<Point*>& pnt_vec)
    : _sfc_pnts(pnt_vec), _bounding_volume(nullptr), _surface_grid(nullptr)
{
}

Surface::Surface(Surface const& src)
    : _sfc_pnts(src._sfc_pnts),
      _bounding_volume(new AABB(*(src._bounding_volume))),
      _surface_grid(nullptr)
{
    _sfc_triangles.reserve(src._sfc_triangles.size());
    std::transform(src._sfc_triangles.cbegin(),
                   src._sfc_triangles.cend(),
                   std::back_inserter(_sfc_triangles),
                   [](Triangle* t) { return new Triangle(*t); });
}

Surface::~Surface()
{
    for (auto& _sfc_triangle : _sfc_triangles)
    {
        delete _sfc_triangle;
    }
}

void Surface::addTriangle(std::size_t pnt_a,
                          std::size_t pnt_b,
                          std::size_t pnt_c)
{
    assert(pnt_a < _sfc_pnts.size() && pnt_b < _sfc_pnts.size() &&
           pnt_c < _sfc_pnts.size());

    // Check if two points of the triangle have identical IDs
    if (pnt_a == pnt_b || pnt_a == pnt_c || pnt_b == pnt_c)
    {
        return;
    }

    // Adding a new triangle invalides the surface grid.
    _surface_grid.reset();

    _sfc_triangles.push_back(new Triangle(_sfc_pnts, pnt_a, pnt_b, pnt_c));
    if (!_bounding_volume)
    {
        std::vector<std::size_t> ids(3);
        ids[0] = pnt_a;
        ids[1] = pnt_b;
        ids[2] = pnt_c;
        _bounding_volume = std::make_unique<GeoLib::AABB>(_sfc_pnts, ids);
    }
    else
    {
        _bounding_volume->update(*_sfc_pnts[pnt_a]);
        _bounding_volume->update(*_sfc_pnts[pnt_b]);
        _bounding_volume->update(*_sfc_pnts[pnt_c]);
    }
}

std::size_t Surface::getNumberOfTriangles() const
{
    return _sfc_triangles.size();
}

const Triangle* Surface::operator[](std::size_t i) const
{
    assert(i < _sfc_triangles.size());
    return _sfc_triangles[i];
}

bool Surface::isPntInBoundingVolume(MathLib::Point3d const& pnt,
                                    double eps) const
{
    return _bounding_volume->containsPoint(pnt, eps);
}

bool Surface::isPntInSfc(MathLib::Point3d const& pnt, double eps) const
{
    // Mutable _surface_grid is constructed if method is called the first time.
    if (_surface_grid == nullptr)
    {
        _surface_grid = std::make_unique<GeoLib::SurfaceGrid>(this);
    }
    return _surface_grid->isPointInSurface(pnt, eps);
}

const Triangle* Surface::findTriangle(MathLib::Point3d const& pnt) const
{
    for (auto _sfc_triangle : _sfc_triangles)
    {
        if (_sfc_triangle->containsPoint(pnt))
        {
            return _sfc_triangle;
        }
    }
    return nullptr;
}

}  // namespace GeoLib
