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

namespace MeshGeoTools
{
class MeshNodesAlongSurface;
}

namespace NumLib
{
/**
 * Class for linearly interpolating values along a GeoLib::Surface object
 */
class LinearInterpolationAlongSurface : public ISpatialFunction
{
public:
	/**
	 * Constructor
	 * @param sfc  a surface object
	 * @param vec_interpolate_point_ids  a vector of point IDs where values are known
	 * @param vec_interpolate_point_values  a vector of values at points
	 * @param default_value  default value when a given point is not located on a surface
	 */
	LinearInterpolationAlongSurface(
			const GeoLib::Surface& sfc,
			const std::vector<std::size_t>& vec_interpolate_point_ids,
			const std::vector<double>& vec_interpolate_point_values,
			const double default_value = 0.0);

	/**
	 * interpolate a value at the given point
	 * @param pnt  a point object
	 * @return interpolated value. A default value is returned if the given point
	 * is not located on a surface
	 */
	double operator()(const GeoLib::Point& pnt) const;

private:
	/// rotate a triangle to XY plane
	void rotate(std::vector<GeoLib::Point> &_pnts) const;

	/// do an interpolation within a triangle
	double interpolateInTri(const GeoLib::Triangle &tri, double const* const vertex_values, double const* const pnt) const;

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

#endif //LINEARINTERPOLATIONSURFACE_H_

