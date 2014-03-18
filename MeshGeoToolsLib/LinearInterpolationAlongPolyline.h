/**
 * \author Norihiro Watanabe
 * \date   2014-03-18
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINEARINTERPOLATIONALONGPOLYLINE_H_
#define LINEARINTERPOLATIONALONGPOLYLINE_H_

#include <vector>

namespace GeoLib
{
class Polyline;
}

namespace MeshGeoTools
{
class MeshNodesAlongPolyline;

/**
 * Class for linearly interpolating values along a GeoLib::Polyline object
 */
class LinearInterpolationAlongPolyline
{
public:
	/**
	 * Linearly interpolate node values from point values along a polyline
	 * @param mshNodesAlongPoly             MeshGeoTools::MeshNodesAlongPolyline object
	 * @param vec_interpolate_point_ids     a vector of point IDs
	 * @param vec_interpolate_point_values  a vector of point values
	 * @param vec_node_values               a vector of interpolated node values
	 */
	LinearInterpolationAlongPolyline(
			const MeshGeoTools::MeshNodesAlongPolyline &mshNodesAlongPoly,
			const std::vector<std::size_t> &vec_interpolate_point_ids,
			const std::vector<double> &vec_interpolate_point_values,
			std::vector<double> &vec_node_values);
};

} // MeshGeoTools

#endif //LINEARINTERPOLATIONALONGPOLYLINE_H_

