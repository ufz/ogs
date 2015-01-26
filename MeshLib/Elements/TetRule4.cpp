/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TetRule4.h"

#include "logog/include/logog.hpp"

#include "GeoLib/AnalyticalGeometry.h"

#include "MeshLib/Node.h"
#include "Tri.h"
#include "Line.h"

namespace MeshLib {

const unsigned TetRule4::n_all_nodes;

const unsigned TetRule4::n_base_nodes;

const unsigned TetRule4::face_nodes[4][3] =
{
	{0, 2, 1}, // Face 0
	{0, 1, 3}, // Face 1
	{1, 2, 3}, // Face 2
	{2, 0, 3}  // Face 3
};

const unsigned TetRule4::edge_nodes[6][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{0, 2}, // Edge 2
	{0, 3}, // Edge 3
	{1, 3}, // Edge 4
	{2, 3}  // Edge 5
};

const Element* TetRule4::getFace(const Element* e, unsigned i)
{
	if (i<e->getNFaces())
	{
		unsigned nFaceNodes (e->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
		return new Tri(nodes);
	}
	ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
	return nullptr;
}

double TetRule4::computeVolume(Node const* const* _nodes)
{
	return GeoLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[3]->getCoords());
}

bool TetRule4::isPntInElement(Node const* const* _nodes, MathLib::Point3d const& pnt, double eps)
{
	return GeoLib::isPointInTetrahedron(pnt, *_nodes[0], *_nodes[1], *_nodes[2], *_nodes[3], eps);
}

unsigned TetRule4::identifyFace(Node const* const* _nodes, Node* nodes[3])
{
	for (unsigned i=0; i<4; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<3; j++)
			for (unsigned k=0; k<3; k++)
				if (_nodes[face_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

ElementErrorCode TetRule4::validate(const Element* e)
{
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();
	error_code[ElementErrorFlag::NodeOrder]  = !e->testElementNodeOrder();
	return error_code;
}

} // end namespace MeshLib
