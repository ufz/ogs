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

#include "LinearInterpolationAlongSurface.h"

#include <limits>
#include <array>
#include <memory>

#include "logog/include/logog.hpp"

#include "MathLib/Vector3.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshGeoToolsLib/MeshNodesAlongSurface.h"

namespace MeshGeoTools
{

LinearInterpolationAlongSurface::LinearInterpolationAlongSurface(
		const MeshLib::Mesh &msh,
		const MeshGeoTools::MeshNodesAlongSurface &mshNodesAlongSurface,
		const std::vector<std::size_t> &vec_interpolate_point_ids,
		const std::vector<double> &vec_interpolate_point_values,
		const double default_value,
		std::vector<double> &vec_node_values)
{
	const std::vector<std::size_t> &nodes_on_sfc = mshNodesAlongSurface.getNodeIDs();
	const GeoLib::Surface &sfc = mshNodesAlongSurface.getSurface();
	vec_node_values.resize(nodes_on_sfc.size());
	std::array<double, 3> pnt_values;
	for (std::size_t i=0; i<vec_node_values.size(); i++ ) {
		const double* node_coords = msh.getNode(nodes_on_sfc[i])->getCoords();
		if (!sfc.isPntInBoundingVolume(node_coords))
			continue;
		auto* tri = sfc.findTriangle(node_coords);
		if (tri == nullptr)
			continue;

		for (unsigned j=0; j<3; j++) {
			auto itr = std::find(vec_interpolate_point_ids.begin(), vec_interpolate_point_ids.end(), (*tri)[j]);
			if (itr != vec_interpolate_point_ids.end()) {
				const std::size_t index = std::distance(vec_interpolate_point_ids.begin(), itr);
				pnt_values[j] = vec_interpolate_point_values[index];
			} else {
				pnt_values[j] = default_value;
			}
		}
		vec_node_values[i] = interpolateInTri(*tri, pnt_values.data(), node_coords);
	}
}

void LinearInterpolationAlongSurface::rotate(std::vector<GeoLib::Point> &pnts) const
{
	std::vector<GeoLib::Point*> p_pnts(3);
	for (unsigned i=0; i<3; i++)
		p_pnts[i] = &pnts[i];

	// calculate supporting plane
	MathLib::Vector3 plane_normal;
	double d;
	// compute the plane normal
	GeoLib::getNewellPlane(p_pnts, plane_normal, d);

	const double tol (std::numeric_limits<double>::epsilon());
	if (std::abs(plane_normal[0]) > tol || std::abs(plane_normal[1]) > tol) {
		// rotate copied points into x-y-plane
		GeoLib::rotatePointsToXY(plane_normal, p_pnts);
	}

	for (std::size_t k(0); k<pnts.size(); k++)
		(*(p_pnts[k]))[2] = 0.0; // should be -= d but there are numerical errors
}


double LinearInterpolationAlongSurface::interpolateInTri(const GeoLib::Triangle &tri, double const* const vertex_values, double const* const pnt) const
{
	std::vector<GeoLib::Point> pnts;
	for (unsigned i=0; i<3; i++)
		pnts.emplace_back(*tri.getPoint(i));
	rotate(pnts);

	const double* v1 = pnts[0].getCoords();
	const double* v2 = pnts[1].getCoords();
	const double* v3 = pnts[2].getCoords();
	const double area = MathLib::calcTriangleArea(v1, v2, v3);

	if (area==.0) {
		double sum = .0;
		for (unsigned i=0; i<3; i++)
			sum += vertex_values[i];
		return sum / static_cast<double>(3);
	}

	double a[3], b[3], c[3];
	// 1st vertex
	a[0] = 0.5/area*(v2[0]*v3[1]-v3[0]*v2[1]);
	b[0] = 0.5/area*(v2[1]-v3[1]);
	c[0] = 0.5/area*(v3[0]-v2[0]);
	// 2nd vertex
	a[1] = 0.5/area*(v3[0]*v1[1]-v1[0]*v3[1]);
	b[1] = 0.5/area*(v3[1]-v1[1]);
	c[1] = 0.5/area*(v1[0]-v3[0]);
	// 3rd vertex
	a[2] = 0.5/area*(v1[0]*v2[1]-v2[0]*v1[1]);
	b[2] = 0.5/area*(v1[1]-v2[1]);
	c[2] = 0.5/area*(v2[0]-v1[0]);

	double val = .0;
	for (unsigned i=0; i<3; i++)
		val += (a[i]+b[i]*pnt[0]+c[i]*pnt[1]) * vertex_values[i];

	return val;
}

} // MeshGeoTools

