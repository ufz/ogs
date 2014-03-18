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

#ifndef LINEARINTERPOLATIONSURFACE_H_
#define LINEARINTERPOLATIONSURFACE_H_

#include <vector>

#include "GeoLib/Point.h"

namespace GeoLib
{
class Surface;
class Triangle;
}

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoTools
{
class MeshNodesAlongSurface;

/**
 * Class for linearly interpolating values along a GeoLib::Surface object
 */
class LinearInterpolationAlongSurface
{
public:
	/**
	 * Linearly interpolate node values from surface point values
	 * @param msh                           a mesh object
	 * @param mshNodesAlongSurface          MeshGeoTools::MeshNodesAlongSurface object
	 * @param vec_interpolate_point_ids     a vector of point IDs
	 * @param vec_interpolate_point_values  a vector of point values
	 * @param default_value                 a default value
	 * @param vec_node_values               a vector of interpolated node values
	 */
	LinearInterpolationAlongSurface(
			const MeshLib::Mesh &msh,
			const MeshGeoTools::MeshNodesAlongSurface &mshNodesAlongSurface,
			const std::vector<std::size_t> &vec_interpolate_point_ids,
			const std::vector<double> &vec_interpolate_point_values,
			const double default_value,
			std::vector<double> &vec_node_values);

private:
	/// rotate a triangle to XY plane
	void rotate(std::vector<GeoLib::Point> &_pnts) const;

	/// do an interpolation within a triangle
	double interpolateInTri(const GeoLib::Triangle &tri, double const* const vertex_values, double const* const pnt) const;

};

} // NumLib

#endif //LINEARINTERPOLATIONSURFACE_H_

