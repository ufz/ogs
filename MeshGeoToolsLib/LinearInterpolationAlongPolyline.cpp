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

#include "LinearInterpolationAlongPolyline.h"

#include "logog/include/logog.hpp"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"

namespace MeshGeoTools
{

LinearInterpolationAlongPolyline::LinearInterpolationAlongPolyline(
		const MeshGeoTools::MeshNodesAlongPolyline &mshNodesAlongPoly,
		const std::vector<std::size_t> &vec_interpolate_point_ids,
		const std::vector<double> &vec_interpolate_point_values,
		std::vector<double> &vec_node_values)
{
	// setup data for interpolation
	std::vector<double> vec_known_dist;
	std::vector<double> vec_known_values;
	vec_known_dist.reserve(vec_interpolate_point_ids.size());
	vec_known_values.reserve(vec_interpolate_point_ids.size());
	const GeoLib::Polyline &ply = mshNodesAlongPoly.getPolyline();
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

	// do interpolation
	vec_node_values.resize(mshNodesAlongPoly.getNodeIDs().size());
	MathLib::PiecewiseLinearInterpolation (
		vec_known_dist,
		vec_known_values,
		mshNodesAlongPoly.getDistOfProjNodeFromPlyStart(),
		vec_node_values);
}

} // MeshGeoTools

