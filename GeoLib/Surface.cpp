/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <list>

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
    : sfc_pnts_(pnt_vec), bounding_volume_(nullptr), surface_grid_(nullptr)
{
}

Surface::Surface(Surface const& src)
    : sfc_pnts_(src.sfc_pnts_),
      bounding_volume_(new AABB(*(src.bounding_volume_))),
      surface_grid_(nullptr)
{
    sfc_triangles_.reserve(src.sfc_triangles_.size());
    std::transform(src.sfc_triangles_.cbegin(),
                   src.sfc_triangles_.cend(),
                   std::back_inserter(sfc_triangles_),
                   [](Triangle* t) { return new Triangle(*t); });
}

Surface::~Surface()
{
    for (auto& sfc_triangle_ : sfc_triangles_)
    {
        delete sfc_triangle_;
    }
}

void Surface::addTriangle(std::size_t pnt_a,
                          std::size_t pnt_b,
                          std::size_t pnt_c)
{
    assert(pnt_a < sfc_pnts_.size() && pnt_b < sfc_pnts_.size() &&
           pnt_c < sfc_pnts_.size());

    // Check if two points of the triangle have identical IDs
    if (pnt_a == pnt_b || pnt_a == pnt_c || pnt_b == pnt_c)
    {
        return;
    }

    // Adding a new triangle invalides the surface grid.
    surface_grid_.reset();

    sfc_triangles_.push_back(new Triangle(sfc_pnts_, pnt_a, pnt_b, pnt_c));
    if (!bounding_volume_)
    {
        std::vector<std::size_t> ids(3);
        ids[0] = pnt_a;
        ids[1] = pnt_b;
        ids[2] = pnt_c;
        bounding_volume_ = std::make_unique<GeoLib::AABB>(sfc_pnts_, ids);
    }
    else
    {
        bounding_volume_->update(*sfc_pnts_[pnt_a]);
        bounding_volume_->update(*sfc_pnts_[pnt_b]);
        bounding_volume_->update(*sfc_pnts_[pnt_c]);
    }
}

std::size_t Surface::getNumberOfTriangles() const
{
    return sfc_triangles_.size();
}

const Triangle* Surface::operator[](std::size_t i) const
{
    assert(i < sfc_triangles_.size());
    return sfc_triangles_[i];
}

bool Surface::isPntInBoundingVolume(MathLib::Point3d const& pnt,
                                    double eps) const
{
    return bounding_volume_->containsPoint(pnt, eps);
}

bool Surface::isPntInSfc(MathLib::Point3d const& pnt, double eps) const
{
    // Mutable surface_grid_ is constructed if method is called the first time.
    if (surface_grid_ == nullptr)
    {
        surface_grid_ = std::make_unique<GeoLib::SurfaceGrid>(this);
    }
    return surface_grid_->isPointInSurface(pnt, eps);
}
}  // namespace GeoLib
