/*
 * \brief Definition of the class MeshElementGrid.
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 */

#include "MeshElementGrid.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <memory>

#include "../Mesh.h"
#include "../Node.h"
#include "../Elements/Element.h"

#include "GeoLib/GEOObjects.h"

namespace MeshLib {

MeshElementGrid::MeshElementGrid(MeshLib::Mesh const& sfc_mesh) :
    aabb_{sfc_mesh.getNodes().cbegin(), sfc_mesh.getNodes().cend()},
    n_steps_({{1,1,1}})
{
    auto getDimensions =
        [](MathLib::Point3d const& min, MathLib::Point3d const& max)
    {
        std::bitset<3> dim;  // all bits are set to zero.
        for (std::size_t k(0); k < 3; ++k) {
            double const tolerance(
                std::nexttoward(max[k],std::numeric_limits<double>::max())-max[k]);
            if (std::abs(max[k] - min[k]) > tolerance)
            {
                dim[k] = true;
            }
        }
        return dim;
    };

    MathLib::Point3d const& min_pnt(aabb_.getMinPoint());
    MathLib::Point3d const& max_pnt(aabb_.getMaxPoint());
    auto const dim = getDimensions(min_pnt, max_pnt);

    std::array<double, 3> delta{{ max_pnt[0] - min_pnt[0],
        max_pnt[1] - min_pnt[1], max_pnt[2] - min_pnt[2] }};

    const std::size_t n_eles(sfc_mesh.getNumberOfElements());
    const std::size_t n_eles_per_cell(100);

    // *** condition: n_eles / n_cells < n_eles_per_cell
    //                where n_cells = n_steps_[0] * n_steps_[1] * n_steps_[2]
    // *** with n_steps_[0] = ceil(pow(n_eles*delta[0]*delta[0]/(n_eles_per_cell*delta[1]*delta[2]), 1/3.)));
    //          n_steps_[1] = n_steps_[0] * delta[1]/delta[0],
    //          n_steps_[2] = n_steps_[0] * delta[2]/delta[0]
    auto sc_ceil = [](double v){
        return static_cast<std::size_t>(std::ceil(v));
    };

    switch (dim.count()) {
    case 3: // 3d case
        n_steps_[0] = sc_ceil(std::cbrt(
            n_eles*delta[0]*delta[0]/(n_eles_per_cell*delta[1]*delta[2])));
        n_steps_[1] = sc_ceil(n_steps_[0] * std::min(delta[1] / delta[0], 100.0));
        n_steps_[2] = sc_ceil(n_steps_[0] * std::min(delta[2] / delta[0], 100.0));
        break;
    case 2: // 2d cases
        if (dim[0] && dim[2]) { // 2d case: xz plane, y = const
            n_steps_[0] = sc_ceil(std::sqrt(n_eles*delta[0]/(n_eles_per_cell*delta[2])));
            n_steps_[2] = sc_ceil(n_steps_[0]*delta[2]/delta[0]);
        }
        else if (dim[0] && dim[1]) { // 2d case: xy plane, z = const
            n_steps_[0] = sc_ceil(std::sqrt(n_eles*delta[0]/(n_eles_per_cell*delta[1])));
            n_steps_[1] = sc_ceil(n_steps_[0] * delta[1]/delta[0]);
        }
        else if (dim[1] && dim[2]) { // 2d case: yz plane, x = const
            n_steps_[1] = sc_ceil(std::sqrt(n_eles*delta[1]/(n_eles_per_cell*delta[2])));
            n_steps_[2] = sc_ceil(n_eles * delta[2] / (n_eles_per_cell*delta[1]));
        }
        break;
    case 1: // 1d cases
        for (std::size_t k(0); k<3; ++k) {
            if (dim[k]) {
                n_steps_[k] = sc_ceil(static_cast<double>(n_eles)/n_eles_per_cell);
            }
        }
    }

    // some frequently used expressions to fill the vector of elements per grid
    // cell
    for (std::size_t k(0); k<3; k++) {
        step_sizes_[k] = delta[k] / n_steps_[k];
        inverse_step_sizes_[k] = 1.0 / step_sizes_[k];
    }

    elements_in_grid_box_.resize(n_steps_[0]*n_steps_[1]*n_steps_[2]);
    sortElementsInGridCells(sfc_mesh);
}

MathLib::Point3d const& MeshElementGrid::getMinPoint() const
{
    return aabb_.getMinPoint();
}

MathLib::Point3d const& MeshElementGrid::getMaxPoint() const
{
    return aabb_.getMaxPoint();
}

void MeshElementGrid::sortElementsInGridCells(MeshLib::Mesh const& sfc_mesh)
{
    for (auto const element : sfc_mesh.getElements()) {
        if (! sortElementInGridCells(*element)) {
            OGS_FATAL("Sorting element (id={:d}) into mesh element grid.",
                      element->getID());
        }
    }
}

bool MeshElementGrid::sortElementInGridCells(MeshLib::Element const& element)
{
    std::array<std::size_t,3> min{};
    std::array<std::size_t,3> max{};
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
        {
            return false;
        }

        for (std::size_t j(0); j < 3; ++j)
        {
            if (min[j] > c.second[j])
            {
                min[j] = c.second[j];
            }
            if (max[j] < c.second[j])
            {
                max[j] = c.second[j];
            }
        }
    }

    const std::size_t n_plane(n_steps_[0]*n_steps_[1]);

    // If a node of an element is almost equal to the upper right point of the
    // AABB the grid cell coordinates computed by getGridCellCoordintes() could
    // be to large (due to numerical errors). The following lines ensure that
    // the grid cell coordinates are in the valid range.
    for (std::size_t k(0); k < 3; ++k)
    {
        max[k] = std::min(n_steps_[k] - 1, max[k]);
    }

    // insert the element into the grid cells
    for (std::size_t i(min[0]); i<=max[0]; i++) {
        for (std::size_t j(min[1]); j<=max[1]; j++) {
            for (std::size_t k(min[2]); k<=max[2]; k++) {
                elements_in_grid_box_[i+j*n_steps_[0]+k*n_plane].push_back(&element);
            }
        }
    }

    return true;
}

std::pair<bool, std::array<std::size_t, 3>>
MeshElementGrid::getGridCellCoordinates(MathLib::Point3d const& p) const
{
    bool valid(true);
    std::array<std::size_t, 3> coords{};

    for (std::size_t k(0); k<3; ++k) {
        const double d(p[k]-aabb_.getMinPoint()[k]);
        if (d < 0.0) {
            valid = false;
            coords[k] = 0;
        } else if (aabb_.getMaxPoint()[k] <= p[k]) {
            valid = false;
            coords[k] = n_steps_[k]-1;
        } else {
            coords[k] = static_cast<std::size_t>(d * inverse_step_sizes_[k]);
        }
    }

    return std::make_pair(valid, coords);
}
}  // namespace MeshLib
