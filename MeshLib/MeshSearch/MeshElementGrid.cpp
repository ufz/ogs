/*
 * \brief Definition of the class MeshElementGrid.
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 */

#include "MeshElementGrid.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <memory>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
MeshElementGrid::MeshElementGrid(MeshLib::Mesh const& mesh)
    : _aabb{mesh.getNodes().cbegin(), mesh.getNodes().cend()},
      _n_steps({{1, 1, 1}})
{
    auto getDimensions = [](auto const& min, auto const& max)
    {
        std::bitset<3> dim;  // all bits are set to zero.
        for (std::size_t k(0); k < 3; ++k)
        {
            double const tolerance(
                std::nexttoward(max[k], std::numeric_limits<double>::max()) -
                max[k]);
            if (std::abs(max[k] - min[k]) > tolerance)
            {
                dim[k] = true;
            }
        }
        return dim;
    };

    auto const& min_pnt(_aabb.getMinPoint());
    auto const& max_pnt(_aabb.getMaxPoint());
    auto const dim = getDimensions(min_pnt, max_pnt);

    std::array<double, 3> delta{{max_pnt[0] - min_pnt[0],
                                 max_pnt[1] - min_pnt[1],
                                 max_pnt[2] - min_pnt[2]}};

    const std::size_t n_eles(mesh.getNumberOfElements());
    const std::size_t n_eles_per_cell(100);

    // *** condition: n_eles / n_cells < n_eles_per_cell
    //                where n_cells = _n_steps[0] * _n_steps[1] * _n_steps[2]
    // *** with _n_steps[0] =
    // ceil(pow(n_eles*delta[0]*delta[0]/(n_eles_per_cell*delta[1]*delta[2]),
    // 1/3.)));
    //          _n_steps[1] = _n_steps[0] * delta[1]/delta[0],
    //          _n_steps[2] = _n_steps[0] * delta[2]/delta[0]
    auto sc_ceil = [](double v)
    { return static_cast<std::size_t>(std::ceil(v)); };

    switch (dim.count())
    {
        case 3:  // 3d case
            _n_steps[0] =
                sc_ceil(std::cbrt(n_eles * delta[0] * delta[0] /
                                  (n_eles_per_cell * delta[1] * delta[2])));
            _n_steps[1] =
                sc_ceil(_n_steps[0] * std::min(delta[1] / delta[0], 100.0));
            _n_steps[2] =
                sc_ceil(_n_steps[0] * std::min(delta[2] / delta[0], 100.0));
            break;
        case 2:  // 2d cases
            if (dim[0] && dim[2])
            {  // 2d case: xz plane, y = const
                _n_steps[0] = sc_ceil(std::sqrt(n_eles * delta[0] /
                                                (n_eles_per_cell * delta[2])));
                _n_steps[2] = sc_ceil(_n_steps[0] * delta[2] / delta[0]);
            }
            else if (dim[0] && dim[1])
            {  // 2d case: xy plane, z = const
                _n_steps[0] = sc_ceil(std::sqrt(n_eles * delta[0] /
                                                (n_eles_per_cell * delta[1])));
                _n_steps[1] = sc_ceil(_n_steps[0] * delta[1] / delta[0]);
            }
            else if (dim[1] && dim[2])
            {  // 2d case: yz plane, x = const
                _n_steps[1] = sc_ceil(std::sqrt(n_eles * delta[1] /
                                                (n_eles_per_cell * delta[2])));
                _n_steps[2] =
                    sc_ceil(n_eles * delta[2] / (n_eles_per_cell * delta[1]));
            }
            break;
        case 1:  // 1d cases
            for (std::size_t k(0); k < 3; ++k)
            {
                if (dim[k])
                {
                    _n_steps[k] =
                        sc_ceil(static_cast<double>(n_eles) / n_eles_per_cell);
                }
            }
    }

    // some frequently used expressions to fill the vector of elements per grid
    // cell
    for (std::size_t k(0); k < 3; k++)
    {
        _step_sizes[k] = delta[k] / _n_steps[k];
        _inverse_step_sizes[k] = 1.0 / _step_sizes[k];
    }

    _elements_in_grid_box.resize(_n_steps[0] * _n_steps[1] * _n_steps[2]);
    sortElementsInGridCells(mesh);
}

Eigen::Vector3d const& MeshElementGrid::getMinPoint() const
{
    return _aabb.getMinPoint();
}

Eigen::Vector3d const& MeshElementGrid::getMaxPoint() const
{
    return _aabb.getMaxPoint();
}

void MeshElementGrid::sortElementsInGridCells(MeshLib::Mesh const& mesh)
{
    for (auto const element : mesh.getElements())
    {
        if (!sortElementInGridCells(*element))
        {
            OGS_FATAL("Sorting element (id={:d}) into mesh element grid.",
                      element->getID());
        }
    }
}

bool MeshElementGrid::sortElementInGridCells(MeshLib::Element const& element)
{
    std::array<std::size_t, 3> min{};
    std::array<std::size_t, 3> max{};
    std::pair<bool, std::array<std::size_t, 3>> c(getGridCellCoordinates(
        *(static_cast<MathLib::Point3d const*>(element.getNode(0)))));
    if (c.first)
    {
        min = c.second;
        max = min;
    }
    else
    {
        return false;
    }

    for (std::size_t k(1); k < element.getNumberOfNodes(); ++k)
    {
        // compute coordinates of the grid for each node of the element
        c = getGridCellCoordinates(
            *(static_cast<MathLib::Point3d const*>(element.getNode(k))));
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

    const std::size_t n_plane(_n_steps[0] * _n_steps[1]);

    // If a node of an element is almost equal to the upper right point of the
    // AABB the grid cell coordinates computed by getGridCellCoordintes() could
    // be to large (due to numerical errors). The following lines ensure that
    // the grid cell coordinates are in the valid range.
    for (std::size_t k(0); k < 3; ++k)
    {
        max[k] = std::min(_n_steps[k] - 1, max[k]);
    }

    // insert the element into the grid cells
    for (std::size_t i(min[0]); i <= max[0]; i++)
    {
        for (std::size_t j(min[1]); j <= max[1]; j++)
        {
            for (std::size_t k(min[2]); k <= max[2]; k++)
            {
                _elements_in_grid_box[i + j * _n_steps[0] + k * n_plane]
                    .push_back(&element);
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

    for (std::size_t k(0); k < 3; ++k)
    {
        const double d(p[k] - _aabb.getMinPoint()[k]);
        if (d < 0.0)
        {
            valid = false;
            coords[k] = 0;
        }
        else if (_aabb.getMaxPoint()[k] <= p[k])
        {
            valid = false;
            coords[k] = _n_steps[k] - 1;
        }
        else
        {
            coords[k] = static_cast<std::size_t>(d * _inverse_step_sizes[k]);
        }
    }

    return std::make_pair(valid, coords);
}
}  // namespace MeshLib
