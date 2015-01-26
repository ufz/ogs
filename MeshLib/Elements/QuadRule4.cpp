/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "QuadRule4.h"

#include "logog/include/logog.hpp"

#include "GeoLib/AnalyticalGeometry.h"

#include "MeshLib/Node.h"

namespace MeshLib {

const unsigned QuadRule4::n_all_nodes;

const unsigned QuadRule4::n_base_nodes;

const unsigned QuadRule4::edge_nodes[4][2] =
{
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{2, 3}, // Edge 2
		{0, 3}  // Edge 3
};

double QuadRule4::computeVolume(Node const* const* _nodes)
{
	return GeoLib::calcTriangleArea(*_nodes[0], *_nodes[1], *_nodes[2])
         + GeoLib::calcTriangleArea(*_nodes[2], *_nodes[3], *_nodes[0]);
}

bool QuadRule4::isPntInElement(Node const* const* _nodes, MathLib::Point3d const& pnt, double eps)
{
	return (GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[1], *_nodes[2], eps) ||
	        GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[2], *_nodes[3], eps));
}

unsigned QuadRule4::identifyFace(Node const* const* _nodes, Node* nodes[3])
{
	for (unsigned i=0; i<4; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<2; j++)
			for (unsigned k=0; k<2; k++)
				if (_nodes[edge_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==2)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

ElementErrorCode QuadRule4::validate(const Element* e)
{
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume]  = e->hasZeroVolume();
	Node const* const* _nodes = e->getNodes();
	error_code[ElementErrorFlag::NonCoplanar] = (!GeoLib::isCoplanar(*_nodes[0], *_nodes[1], *_nodes[2], *_nodes[3]));
	// for collapsed quads (i.e. reduced to a line) this test might result "false" as all four points are actually located on a line.
	if (!error_code[ElementErrorFlag::ZeroVolume])
		error_code[ElementErrorFlag::NonConvex]   = (!(GeoLib::dividedByPlane(*_nodes[0], *_nodes[2], *_nodes[1], *_nodes[3]) &&
			                                           GeoLib::dividedByPlane(*_nodes[1], *_nodes[3], *_nodes[0], *_nodes[2])));
	error_code[ElementErrorFlag::NodeOrder]  = !e->testElementNodeOrder();
	return error_code;
}

} // end namespace MeshLib
