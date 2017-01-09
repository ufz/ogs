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

#include "LinearInterpolationAlongPolyline.h"

#include "GeoLib/Polyline.h"
#include "MeshLib/Mesh.h"
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"

namespace NumLib
{

LinearInterpolationAlongPolyline::LinearInterpolationAlongPolyline(
        const GeoLib::Polyline& ply,
        const std::vector<std::size_t>& vec_interpolate_point_ids,
        const std::vector<double>& vec_interpolate_point_values,
        const double search_length, const double default_value)
: _ply(ply),
  _interpolation(createInterpolation(ply, vec_interpolate_point_ids, vec_interpolate_point_values)),
 _search_length(search_length), _default_value(default_value)
{}

MathLib::PiecewiseLinearInterpolation LinearInterpolationAlongPolyline::createInterpolation(
        const GeoLib::Polyline& ply,
        const std::vector<std::size_t>& vec_interpolate_point_ids,
        const std::vector<double>& vec_interpolate_point_values)
{
    std::vector<double> vec_known_dist;
    std::vector<double> vec_known_values;
    vec_known_dist.reserve(vec_interpolate_point_ids.size());
    vec_known_values.reserve(vec_interpolate_point_ids.size());
    for (std::size_t i=0; i<vec_interpolate_point_ids.size(); i++)
    {
        const std::size_t pnt_id = vec_interpolate_point_ids[i];
        if (!ply.isPointIDInPolyline(pnt_id))
            continue;

        for (std::size_t j=0; j<ply.getNumberOfPoints(); j++)
        {
            if (pnt_id == ply.getPointID(j))
            {
                vec_known_dist.push_back(ply.getLength(j));
                vec_known_values.push_back(vec_interpolate_point_values[i]);
                break;
            }
        }
    }

    return MathLib::PiecewiseLinearInterpolation{std::move(vec_known_dist),
                                                 std::move(vec_known_values)};
}

double LinearInterpolationAlongPolyline::operator()(const MathLib::Point3d& pnt) const
{
    const double dist = _ply.getDistanceAlongPolyline(pnt, _search_length);
    return dist>=0 ? _interpolation.getValue(dist) : _default_value;
}

} // NumLib

