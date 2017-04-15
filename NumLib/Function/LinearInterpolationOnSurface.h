/**
 * \author Norihiro Watanabe
 * \date   2014-03-18
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "GeoLib/Point.h"
#include "ISpatialFunction.h"

namespace GeoLib
{
class Surface;
class Triangle;
}

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
/**
 * Class for linearly interpolating values along a GeoLib::Surface object
 */
class LinearInterpolationOnSurface : public ISpatialFunction
{
public:
    /**
     * Constructor
     * @param sfc  a surface object
     * @param vec_interpolate_point_ids  a vector of point IDs where values are known
     * @param vec_interpolate_point_values  a vector of values at points
     * @param default_value  default value when a given point is not located on a surface
     */
    LinearInterpolationOnSurface(
            const GeoLib::Surface& sfc,
            const std::vector<std::size_t>& vec_interpolate_point_ids,
            const std::vector<double>& vec_interpolate_point_values,
            const double default_value);

    /**
     * interpolate a value at the given point
     * @param pnt  a point object
     * @return interpolated value. A default value is returned if the given point
     * is not located on a surface
     */
    double operator()(const MathLib::Point3d& pnt) const override;

private:
    /// rotate a triangle to XY plane
    void rotate(std::vector<GeoLib::Point> &_pnts) const;

    /**
     * do an interpolation within a triangle. Interpolation is done by shape functions
     * for triangle elements based on areal coordinates as
     * \f[
     *   u(x,y) = sum_i[N_i(x,y)*u_i]  (i=1,2,3)
     * \f]
     * where \f$N_i\f$ is a shape function for node i, \f$u_i\f$ is a value at node i.
     * The shape function is given by
     * \f[
     *   N_i(x,y) = 1/(2A)*(a_i + b_i*x + c_i*y)
     * \f]
     * with an element area \f$A\f$ and geometric parameters \f$a_i\f$, \f$b_i\f$ and \f$c_i\f$.
     * reference: http://kratos-wiki.cimne.upc.edu/index.php/Two-dimensional_Shape_Functions
     *
     * @param tri
     * @param vertex_values
     * @param pnt
     * @return
     */
    double interpolateInTri(const GeoLib::Triangle &tri,
        double const* const vertex_values,
        MathLib::Point3d const& pnt) const;

    /// a surface object
    const GeoLib::Surface& _sfc;
    /// a vector of point IDs where values are known
    const std::vector<std::size_t>& _vec_interpolate_point_ids;
    /// a vector of values at points
    const std::vector<double>& _vec_interpolate_point_values;
    /// default value
    const double _default_value;
};

} // NumLib
