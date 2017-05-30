/*
 * \brief Definition of the class MeshElementGrid.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 */

#include "MeshElementGrid.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <memory>

#include <logog/include/logog.hpp>

#include "../Mesh.h"
#include "../Node.h"
#include "../Elements/Element.h"

#include "GeoLib/GEOObjects.h"

namespace MeshLib {

MeshElementGrid::MeshElementGrid(MeshLib::Mesh const& sfc_mesh) :
    _aabb{sfc_mesh.getNodes().cbegin(), sfc_mesh.getNodes().cend()},
    _n_steps({{1,1,1}})
{
    auto getDimensions =
        [](MathLib::Point3d const& min, MathLib::Point3d const& max)
    {
        std::bitset<3> dim;  // all bits are set to zero.
        for (std::size_t k(0); k < 3; ++k) {
            double const tolerance(
                std::nexttoward(max[k],std::numeric_limits<double>::max())-max[k]);
            if (std::abs(max[k]-min[k]) > tolerance)
                dim[k] = true;
        }
        return dim;
    };

    MathLib::Point3d const& min_pnt(_aabb.getMinPoint());
    MathLib::Point3d const& max_pnt(_aabb.getMaxPoint());
    auto const dim = getDimensions(min_pnt, max_pnt);

    std::array<double, 3> delta{{ max_pnt[0] - min_pnt[0],
        max_pnt[1] - min_pnt[1], max_pnt[2] - min_pnt[2] }};

    const std::size_t n_eles(sfc_mesh.getNumberOfElements());
    const std::size_t n_eles_per_cell(100);

    // *** condition: n_eles / n_cells < n_eles_per_cell
    //                where n_cells = _n_steps[0] * _n_steps[1] * _n_steps[2]
    // *** with _n_steps[0] = ceil(pow(n_eles*delta[0]*delta[0]/(n_eles_per_cell*delta[1]*delta[2]), 1/3.)));
    //          _n_steps[1] = _n_steps[0] * delta[1]/delta[0],
    //          _n_steps[2] = _n_steps[0] * delta[2]/delta[0]
    auto sc_ceil = [](double v){
        return static_cast<std::size_t>(std::ceil(v));
    };

    switch (dim.count()) {
    case 3: // 3d case
        _n_steps[0] = sc_ceil(std::cbrt(
            n_eles*delta[0]*delta[0]/(n_eles_per_cell*delta[1]*delta[2])));
        _n_steps[1] = sc_ceil(_n_steps[0] * std::min(delta[1] / delta[0], 100.0));
        _n_steps[2] = sc_ceil(_n_steps[0] * std::min(delta[2] / delta[0], 100.0));
        break;
    case 2: // 2d cases
        if (dim[0] && dim[2]) { // 2d case: xz plane, y = const
            _n_steps[0] = sc_ceil(std::sqrt(n_eles*delta[0]/(n_eles_per_cell*delta[2])));
            _n_steps[2] = sc_ceil(_n_steps[0]*delta[2]/delta[0]);
        }
        else if (dim[0] && dim[1]) { // 2d case: xy plane, z = const
            _n_steps[0] = sc_ceil(std::sqrt(n_eles*delta[0]/(n_eles_per_cell*delta[1])));
            _n_steps[1] = sc_ceil(_n_steps[0] * delta[1]/delta[0]);
        }
        else if (dim[1] && dim[2]) { // 2d case: yz plane, x = const
            _n_steps[1] = sc_ceil(std::sqrt(n_eles*delta[1]/(n_eles_per_cell*delta[2])));
            _n_steps[2] = sc_ceil(n_eles * delta[2] / (n_eles_per_cell*delta[1]));
        }
        break;
    case 1: // 1d cases
        for (std::size_t k(0); k<3; ++k) {
            if (dim[k]) {
                _n_steps[k] = sc_ceil(static_cast<double>(n_eles)/n_eles_per_cell);
            }
        }
    }

    // some frequently used expressions to fill the vector of elements per grid
    // cell
    for (std::size_t k(0); k<3; k++) {
        _step_sizes[k] = delta[k] / _n_steps[k];
        _inverse_step_sizes[k] = 1.0 / _step_sizes[k];
    }

    _elements_in_grid_box.resize(_n_steps[0]*_n_steps[1]*_n_steps[2]);
    sortElementsInGridCells(sfc_mesh);
}

MathLib::Point3d const& MeshElementGrid::getMinPoint() const
{
    return _aabb.getMinPoint();
}

MathLib::Point3d const& MeshElementGrid::getMaxPoint() const
{
    return _aabb.getMaxPoint();
}

void MeshElementGrid::sortElementsInGridCells(MeshLib::Mesh const& sfc_mesh)
{
    for (auto const element : sfc_mesh.getElements()) {
        if (! sortElementInGridCells(*element)) {
            OGS_FATAL("Sorting element (id=%d) into mesh element grid.",
                element->getID());
        }
    }
}

bool MeshElementGrid::sortElementInGridCells(MeshLib::Element const& element)
{
    std::array<std::size_t,3> min;
    std::array<std::size_t,3> max;
    std::pair<bool, std::array<std::size_t, 3>> c(
        getGridCellCoordinates(*(static_cast<MathLib::Point3d const*>(element.getNode(0)))));
    if (c.first) {
        min = c.second;
        max = min;
    } else {
        return false;
    }

    std::vector<std::array<std::size_t,3>> coord_vecs(element.getNumberOfNodes());
    for (std::size_t k(1); k<element.getNumberOfNodes(); ++k) {
        // compute coordinates of the grid for each node of the element
        c = getGridCellCoordinates(*(static_cast<MathLib::Point3d const*>(element.getNode(k))));
        if (!c.first)
            return false;

        for (std::size_t j(0); j < 3; ++j)
        {
            if (min[j] > c.second[j])
                min[j] = c.second[j];
            if (max[j] < c.second[j])
                max[j] = c.second[j];
        }
    }

    const std::size_t n_plane(_n_steps[0]*_n_steps[1]);

    // If a node of an element is almost equal to the upper right point of the
    // AABB the grid cell coordinates computed by getGridCellCoordintes() could
    // be to large (due to numerical errors). The following lines ensure that
    // the grid cell coordinates are in the valid range.
    for (std::size_t k(0); k<3; ++k)
        max[k] = std::min(_n_steps[k]-1, max[k]);

    // insert the element into the grid cells
    for (std::size_t i(min[0]); i<=max[0]; i++) {
        for (std::size_t j(min[1]); j<=max[1]; j++) {
            for (std::size_t k(min[2]); k<=max[2]; k++) {
                _elements_in_grid_box[i+j*_n_steps[0]+k*n_plane].push_back(&element);
            }
        }
    }

    return true;
}

std::pair<bool, std::array<std::size_t, 3>>
MeshElementGrid::getGridCellCoordinates(MathLib::Point3d const& p) const
{
    bool valid(true);
    std::array<std::size_t, 3> coords;

    for (std::size_t k(0); k<3; ++k) {
        const double d(p[k]-_aabb.getMinPoint()[k]);
        if (d < 0.0) {
            valid = false;
            coords[k] = 0;
        } else if (_aabb.getMaxPoint()[k] <= p[k]) {
            valid = false;
            coords[k] = _n_steps[k]-1;
        } else {
            coords[k] = static_cast<std::size_t>(d * _inverse_step_sizes[k]);
        }
    }

    return std::make_pair(valid, coords);
}

#ifndef NDEBUG
void getGridGeometry(MeshElementGrid const& grid,
                     GeoLib::GEOObjects& geometries,
                     std::string& geometry_name)
{
    std::vector<std::string> cell_names;

    auto addPoints = [&geometries](
        MathLib::Point3d const& p, std::array<double, 3> const& d,
        std::array<std::size_t, 3> const& c, std::string& name) {
        auto pnts = std::make_unique<std::vector<GeoLib::Point*>>();
        pnts->push_back(new GeoLib::Point(p[0]+c[0]*d[0], p[1]+c[1]*d[1], p[2]+c[2]*d[2]));
        pnts->push_back(new GeoLib::Point(p[0]+c[0]*d[0], p[1]+(c[1]+1)*d[1], p[2]+c[2]*d[2]));
        pnts->push_back(new GeoLib::Point(p[0]+(c[0]+1)*d[0], p[1]+(c[1]+1)*d[1], p[2]+c[2]*d[2]));
        pnts->push_back(new GeoLib::Point(p[0]+(c[0]+1)*d[0], p[1]+c[1]*d[1], p[2]+c[2]*d[2]));
        pnts->push_back(new GeoLib::Point(p[0]+c[0]*d[0], p[1]+c[1]*d[1], p[2]+(c[2]+1)*d[2]));
        pnts->push_back(new GeoLib::Point(p[0]+c[0]*d[0], p[1]+(c[1]+1)*d[1], p[2]+(c[2]+1)*d[2]));
        pnts->push_back(new GeoLib::Point(p[0]+(c[0]+1)*d[0], p[1]+(c[1]+1)*d[1], p[2]+(c[2]+1)*d[2]));
        pnts->push_back(new GeoLib::Point(p[0]+(c[0]+1)*d[0], p[1]+c[1]*d[1], p[2]+(c[2]+1)*d[2]));
        std::array<double,3> ulps; // unit in the last place
        double const towards(std::numeric_limits<double>::max());
        ulps[0] = std::nextafter(d[0], towards)-d[0];
        ulps[1] = std::nextafter(d[1], towards)-d[1];
        ulps[2] = std::nextafter(d[2], towards)-d[2];
        double const tolerance(std::min(std::min(ulps[0], ulps[1]), ulps[2]));
        geometries.addPointVec(std::move(pnts), name, nullptr, tolerance);
    };

    for (std::size_t i(0); i<grid._n_steps[0]; ++i) {
        for (std::size_t j(0); j<grid._n_steps[1]; ++j) {
            for (std::size_t k(0); k<grid._n_steps[2]; ++k) {
                cell_names.emplace_back("Grid-"+std::to_string(i)+"-"
                    +std::to_string(j)+"-"+std::to_string(k));
                addPoints(grid._aabb.getMinPoint(), grid._step_sizes,
                          {{i, j, k}}, cell_names.back());
                auto plys = std::make_unique<std::vector<GeoLib::Polyline*>>();
                auto & points = *geometries.getPointVec(cell_names.back());

                auto* ply_bottom(new GeoLib::Polyline(points));
                for (std::size_t l(0); l < 4; ++l)
                    ply_bottom->addPoint(l);
                ply_bottom->addPoint(0); // close to bottom surface
                plys->push_back(ply_bottom);

                auto* ply_top(new GeoLib::Polyline(points));
                for (std::size_t l(4); l<8; ++l)
                    ply_top->addPoint(l);
                ply_top->addPoint(4); // close to top surface
                plys->push_back(ply_top);

                auto* ply_04(new GeoLib::Polyline(points));
                ply_04->addPoint(0);
                ply_04->addPoint(4);
                plys->push_back(ply_04);

                auto* ply_15(new GeoLib::Polyline(points));
                ply_15->addPoint(1);
                ply_15->addPoint(5);
                plys->push_back(ply_15);

                auto* ply_26(new GeoLib::Polyline(points));
                ply_26->addPoint(2);
                ply_26->addPoint(6);
                plys->push_back(ply_26);

                auto* ply_37(new GeoLib::Polyline(points));
                ply_37->addPoint(3);
                ply_37->addPoint(7);
                plys->push_back(ply_37);

                geometries.addPolylineVec(std::move(plys), cell_names.back(),
                    nullptr);
            }
        }
    }
    if (geometries.mergeGeometries(cell_names, geometry_name) == 2)
        geometry_name = cell_names.front();
}
#endif // NDEBUG

} // end namespace MeshLib
