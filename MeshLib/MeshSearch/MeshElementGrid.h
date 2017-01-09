/*
 * \brief Declaration of the MeshElementGrid class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
#include <limits>
#include <vector>

#include "GeoLib/AABB.h"
#include "MathLib/Point3d.h"

namespace GeoLib {
class GEOObjects;
}

namespace MeshLib {

// forward declarations
class Mesh;
class Element;

/// MeshElementGrid implements a grid data structure supporting search
/// operations and covers a given mesh. It consists of grid cells that are all
/// of the same size. Grid cells contain pointers to intersecting mesh elements.
/// @attention The user has to ensure the validity of the pointers while the
/// MeshElementGrid instance lives.
class MeshElementGrid final {
#ifndef NDEBUG
    friend void getGridGeometry(MeshLib::MeshElementGrid const& grid,
                                GeoLib::GEOObjects& geometries,
                                std::string& geometry_name);
#endif
public:
    /// Constructs a grid. Grid cells contains intersecting mesh elements.
    /// @param mesh the MeshLib::Mesh instance the grid will be constructed from
    explicit MeshElementGrid(MeshLib::Mesh const& mesh);

    /// Fill and return a vector containing elements of all grid cells that have
    /// a non-empty intersection with the box that is defined by min and max.
    /// @param min min point of the box
    /// @param max max point of the box
    /// @return a (possible empty) vector of elements
    template <typename POINT>
    std::vector<MeshLib::Element const*> getElementsInVolume(
        POINT const& min, POINT const& max) const
    {
        auto const min_coords(getGridCellCoordinates(min));
        auto const max_coords(getGridCellCoordinates(max));

        std::vector<MeshLib::Element const*> elements_vec;

        const std::size_t n_plane(_n_steps[0]*_n_steps[1]);
        for (std::size_t i(min_coords.second[0]); i<=max_coords.second[0]; i++) {
            for (std::size_t j(min_coords.second[1]); j<=max_coords.second[1]; j++) {
                for (std::size_t k(min_coords.second[2]); k<=max_coords.second[2]; k++) {
                    std::size_t idx(i+j*_n_steps[0]+k*n_plane);
                    std::copy(_elements_in_grid_box[idx].begin(),
                              _elements_in_grid_box[idx].end(),
                              std::back_inserter(elements_vec));
                }
            }
        }
        return elements_vec;
    }

    /// Returns the min point of the internal AABB. The method is a wrapper for
    /// GeoLib::AABB::getMinPoint().
    MathLib::Point3d const& getMinPoint() const;
    /// Returns the max point of the internal AABB. The method is a wrapper for
    /// AABB::getMaxPoint().
    MathLib::Point3d const& getMaxPoint() const;

private:
    void sortElementsInGridCells(MeshLib::Mesh const& sfc_mesh);
    bool sortElementInGridCells(MeshLib::Element const& element);

    GeoLib::AABB _aabb;
    /// Computes the grid cell coordinates for given point. The first element of
    /// the returned pair (bool) is true if the point is within the grid, else
    /// false.
    std::pair<bool, std::array<std::size_t,3>>
        getGridCellCoordinates(MathLib::Point3d const& p) const;
    std::array<double,3> _step_sizes;
    std::array<double,3> _inverse_step_sizes;
    std::array<std::size_t,3> _n_steps;
    std::vector<std::vector<MeshLib::Element const*>> _elements_in_grid_box;
};

#ifndef NDEBUG
/// Transfers the grid cells to a geometry. The grid-geometry can be viewed
/// using the DataExplorer. The visualization can be used to explain the
/// algorithm or to find bugs in the algorithm.
void getGridGeometry(MeshLib::MeshElementGrid const& grid,
                     GeoLib::GEOObjects& geometries,
                     std::string& geometry_name);
#endif

} // end namespace MeshLib
