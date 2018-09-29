/*
 * \date 2012-09-22
 * \brief Definition of the SurfaceGrid class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 */

#include "SurfaceGrid.h"

#include <algorithm>
#include <bitset>

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"

#include "MathLib/Point3d.h"

#include "Surface.h"
#include "Triangle.h"

namespace GeoLib {

SurfaceGrid::SurfaceGrid(Surface const*const sfc) :
    AABB(sfc->getAABB()), _n_steps({{1,1,1}})
{
    // enlarge the bounding, such that the points with maximal coordinates
    // fits into the grid
    for (std::size_t k(0); k<3; ++k) {
        _max_pnt[k] += std::abs(_max_pnt[k]) * 1e-6;
        if (std::abs(_max_pnt[k]) < std::numeric_limits<double>::epsilon()) {
            _max_pnt[k] = (_max_pnt[k] - _min_pnt[k]) * (1.0 + 1e-6);
        }
    }

    std::array<double, 3> delta{{ _max_pnt[0] - _min_pnt[0],
        _max_pnt[1] - _min_pnt[1], _max_pnt[2] - _min_pnt[2] }};

    if (delta[0] < std::numeric_limits<double>::epsilon()) {
        const double max_delta(std::max(delta[1], delta[2]));
        _min_pnt[0] -= max_delta * 0.5e-3;
        _max_pnt[0] += max_delta * 0.5e-3;
        delta[0] = _max_pnt[0] - _min_pnt[0];
    }

    if (delta[1] < std::numeric_limits<double>::epsilon()) {
        const double max_delta(std::max(delta[0], delta[2]));
        _min_pnt[1] -= max_delta * 0.5e-3;
        _max_pnt[1] += max_delta * 0.5e-3;
        delta[1] = _max_pnt[1] - _min_pnt[1];
    }

    if (delta[2] < std::numeric_limits<double>::epsilon()) {
        const double max_delta(std::max(delta[0], delta[1]));
        _min_pnt[2] -= max_delta * 0.5e-3;
        _max_pnt[2] += max_delta * 0.5e-3;
        delta[2] = _max_pnt[2] - _min_pnt[2];
    }

    const std::size_t n_tris(sfc->getNumberOfTriangles());
    const std::size_t n_tris_per_cell(5);

    std::bitset<3> dim; // all bits are set to zero.
    for (std::size_t k(0); k<3; ++k)
        if (std::abs(delta[k]) >= std::numeric_limits<double>::epsilon())
            dim[k] = true;

    // *** condition: n_tris / n_cells < n_tris_per_cell
    //                where n_cells = _n_steps[0] * _n_steps[1] * _n_steps[2]
    // *** with _n_steps[0] = ceil(pow(n_tris*delta[0]*delta[0]/(n_tris_per_cell*delta[1]*delta[2]), 1/3.)));
    //          _n_steps[1] = _n_steps[0] * delta[1]/delta[0],
    //          _n_steps[2] = _n_steps[0] * delta[2]/delta[0]
    auto sc_ceil = [](double v){
        return static_cast<std::size_t>(std::ceil(v));
    };
    switch (dim.count()) {
    case 3: // 3d case
        _n_steps[0] = sc_ceil(std::cbrt(
            n_tris*delta[0]*delta[0]/(n_tris_per_cell*delta[1]*delta[2])));
        _n_steps[1] = sc_ceil(_n_steps[0] * std::min(delta[1] / delta[0], 100.0));
        _n_steps[2] = sc_ceil(_n_steps[0] * std::min(delta[2] / delta[0], 100.0));
        break;
    case 2: // 2d cases
        if (dim[0] && dim[2]) { // 2d case: xz plane, y = const
            _n_steps[0] = sc_ceil(std::sqrt(n_tris*delta[0]/(n_tris_per_cell*delta[2])));
            _n_steps[2] = sc_ceil(_n_steps[0]*delta[2]/delta[0]);
        }
        else if (dim[0] && dim[1]) { // 2d case: xy plane, z = const
            _n_steps[0] = sc_ceil(std::sqrt(n_tris*delta[0]/(n_tris_per_cell*delta[1])));
            _n_steps[1] = sc_ceil(_n_steps[0] * delta[1]/delta[0]);
        }
        else if (dim[1] && dim[2]) { // 2d case: yz plane, x = const
            _n_steps[1] = sc_ceil(std::sqrt(n_tris*delta[1]/(n_tris_per_cell*delta[2])));
            _n_steps[2] = sc_ceil(n_tris * delta[2] / (n_tris_per_cell*delta[1]));
        }
        break;
    case 1: // 1d cases
        for (std::size_t k(0); k<3; ++k) {
            if (dim[k]) {
                _n_steps[k] = sc_ceil(static_cast<double>(n_tris)/n_tris_per_cell);
            }
        }
    }

    // some frequently used expressions to fill the grid vectors
    for (std::size_t k(0); k<3; k++) {
        _step_sizes[k] = delta[k] / _n_steps[k];
        _inverse_step_sizes[k] = 1.0 / _step_sizes[k];
    }

    _triangles_in_grid_box.resize(_n_steps[0]*_n_steps[1]*_n_steps[2]);
    sortTrianglesInGridCells(sfc);
}

void SurfaceGrid::sortTrianglesInGridCells(Surface const*const sfc)
{
    for (std::size_t l(0); l<sfc->getNumberOfTriangles(); l++) {
        if (! sortTriangleInGridCells((*sfc)[l])) {
            Point const& p0(*((*sfc)[l]->getPoint(0)));
            Point const& p1(*((*sfc)[l]->getPoint(1)));
            Point const& p2(*((*sfc)[l]->getPoint(2)));
            ERR("Sorting triangle %d [(%f,%f,%f), (%f,%f,%f), (%f,%f,%f) into "
                "grid.",
                l, p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]
            );
            OGS_FATAL("");
        }
    }
}

bool SurfaceGrid::sortTriangleInGridCells(Triangle const*const triangle)
{
    // compute grid coordinates for each triangle point
    boost::optional<std::array<std::size_t, 3> const> c_p0(
        getGridCellCoordinates(*(triangle->getPoint(0))));
    if (!c_p0)
        return false;
    boost::optional<std::array<std::size_t, 3> const> c_p1(
        getGridCellCoordinates(*(triangle->getPoint(1))));
    if (!c_p1)
        return false;
    boost::optional<std::array<std::size_t, 3> const> c_p2(
        getGridCellCoordinates(*(triangle->getPoint(2))));
    if (!c_p2)
        return false;

    // determine interval in grid (grid cells the triangle will be inserted)
    std::size_t const i_min(std::min(std::min((*c_p0)[0], (*c_p1)[0]), (*c_p2)[0]));
    std::size_t const i_max(std::max(std::max((*c_p0)[0], (*c_p1)[0]), (*c_p2)[0]));
    std::size_t const j_min(std::min(std::min((*c_p0)[1], (*c_p1)[1]), (*c_p2)[1]));
    std::size_t const j_max(std::max(std::max((*c_p0)[1], (*c_p1)[1]), (*c_p2)[1]));
    std::size_t const k_min(std::min(std::min((*c_p0)[2], (*c_p1)[2]), (*c_p2)[2]));
    std::size_t const k_max(std::max(std::max((*c_p0)[2], (*c_p1)[2]), (*c_p2)[2]));

    const std::size_t n_plane(_n_steps[0]*_n_steps[1]);

    // insert the triangle into the grid cells
    for (std::size_t i(i_min); i<=i_max; i++) {
        for (std::size_t j(j_min); j<=j_max; j++) {
            for (std::size_t k(k_min); k<=k_max; k++) {
                _triangles_in_grid_box[i+j*_n_steps[0]+k*n_plane].push_back(triangle);
            }
        }
    }

    return true;
}

boost::optional<std::array<std::size_t, 3>>
SurfaceGrid::getGridCellCoordinates(MathLib::Point3d const& p) const
{
    std::array<std::size_t, 3> coords{{
        static_cast<std::size_t>((p[0]-_min_pnt[0]) * _inverse_step_sizes[0]),
        static_cast<std::size_t>((p[1]-_min_pnt[1]) * _inverse_step_sizes[1]),
        static_cast<std::size_t>((p[2]-_min_pnt[2]) * _inverse_step_sizes[2])
    }};

    if (coords[0]>=_n_steps[0] || coords[1]>=_n_steps[1] || coords[2]>=_n_steps[2]) {
        DBUG("Computed indices (%d,%d,%d), max grid cell indices (%d,%d,%d)",
            coords[0], coords[1], coords[2], _n_steps[0], _n_steps[1], _n_steps[2]);
        return boost::optional<std::array<std::size_t, 3>>();
    }
    return boost::optional<std::array<std::size_t, 3>>(coords);
}

bool SurfaceGrid::isPointInSurface(MathLib::Point3d const& p, double eps) const
{
    boost::optional<std::array<std::size_t,3>> optional_c(getGridCellCoordinates(p));
    if (!optional_c) {
        return false;
    }
    std::array<std::size_t,3> c(optional_c.get());

    std::size_t const grid_cell_idx(c[0]+c[1]*_n_steps[0]+c[2]*_n_steps[0]*_n_steps[1]);
    std::vector<Triangle const*> const& triangles(_triangles_in_grid_box[grid_cell_idx]);
    bool nfound(true);
    const std::size_t n_triangles(triangles.size());
    for (std::size_t k(0); k < n_triangles && nfound; k++) {
        if (triangles[k]->containsPoint(p, eps)) {
            nfound = false;
        }
    }

    return !nfound;
}

} // end namespace GeoLib
