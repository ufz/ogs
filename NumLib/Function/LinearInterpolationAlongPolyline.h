/**
 * \author Norihiro Watanabe
 * \date   2014-03-18
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "GeoLib/Point.h"
#include "ISpatialFunction.h"

namespace GeoLib
{
class Polyline;
}

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
/**
 * Class for linearly interpolating values along a GeoLib::Polyline object
 */
class LinearInterpolationAlongPolyline : public ISpatialFunction
{
public:
    /**
     * Constructor
     * @param ply  a polyline object
     * @param vec_interpolate_point_ids  a vector of point IDs where values are known
     * @param vec_interpolate_point_values  a vector of values at points
     * @param search_length  distance threshold to decide whether a point is located on a polyline or not
     * @param default_value  default value when a given point is not located on a polyine
     */
    LinearInterpolationAlongPolyline(
            const GeoLib::Polyline& ply,
            const std::vector<std::size_t>& vec_interpolate_point_ids,
            const std::vector<double>& vec_interpolate_point_values,
            const double search_length, const double default_value);

    /**
     * interpolate a value at the given point
     * @param pnt  a point object
     * @return interpolated value. A default value is returned if the given point
     * is not located on a polyline
     */
    double operator()(const MathLib::Point3d& pnt) const override;

private:
    /// construct an interpolation algorithm
    static MathLib::PiecewiseLinearInterpolation createInterpolation(
            const GeoLib::Polyline& ply,
            const std::vector<std::size_t>& vec_interpolate_point_ids,
            const std::vector<double>& vec_interpolate_point_values);

    /// a polyline object
    const GeoLib::Polyline& _ply;
    /// an interpolation algorithm
    const MathLib::PiecewiseLinearInterpolation _interpolation;
    /// search length
    const double _search_length;
    /// default value
    const double _default_value;
};

} // NumLib
