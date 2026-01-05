// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>
#include <limits>
#include <optional>
#include <vector>

#include "AABB.h"
#include "Point.h"

namespace GeoLib
{

// forward declarations
class Triangle;
class Surface;

class SurfaceGrid final : public AABB
{
public:
    explicit SurfaceGrid(GeoLib::Surface const* const sfc);
    bool isPointInSurface(
        MathLib::Point3d const& pnt,
        double eps = std::numeric_limits<double>::epsilon()) const;

private:
    void sortTrianglesInGridCells(GeoLib::Surface const* const sfc);
    bool sortTriangleInGridCells(GeoLib::Triangle const* const triangle);
    std::optional<std::array<std::size_t, 3>> getGridCellCoordinates(
        MathLib::Point3d const& p) const;
    std::array<double, 3> _step_sizes{};
    std::array<double, 3> _inverse_step_sizes{};
    std::array<std::size_t, 3> _n_steps;
    std::vector<std::vector<GeoLib::Triangle const*>> _triangles_in_grid_box;
};

}  // end namespace GeoLib
