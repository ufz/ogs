/*
 * \date 2012-09-22
 * \brief Declaration of the SurfaceGrid class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef SURFACEGRID_H_
#define SURFACEGRID_H_

#include <array>
#include <limits>
#include <vector>

#include <boost/optional.hpp>

#include "AABB.h"
#include "Point.h"

namespace GeoLib {

// forward declarations
class Triangle;
class Surface;

class SurfaceGrid : public AABB {
public:
    explicit SurfaceGrid(GeoLib::Surface const*const sfc);
    bool isPointInSurface(MathLib::Point3d const & pnt,
        double eps = std::numeric_limits<double>::epsilon()) const;

private:
    void sortTrianglesInGridCells(GeoLib::Surface const*const surface);
    bool sortTriangleInGridCells(GeoLib::Triangle const*const triangle);
    boost::optional<std::array<std::size_t,3>>
        getGridCellCoordinates(MathLib::Point3d const& p) const;
    std::array<double,3> _step_sizes;
    std::array<double,3> _inverse_step_sizes;
    std::array<std::size_t,3> _n_steps;
    std::vector<std::vector<GeoLib::Triangle const*>> _triangles_in_grid_box;
};

} // end namespace GeoLib

#endif /* SURFACEGRID_H_ */
