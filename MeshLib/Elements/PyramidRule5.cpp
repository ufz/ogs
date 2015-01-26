/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PyramidRule5.h"

#include "logog/include/logog.hpp"

#include "GeoLib/AnalyticalGeometry.h"

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib {

const unsigned PyramidRule5::n_all_nodes;

const unsigned PyramidRule5::n_base_nodes;

const unsigned PyramidRule5::face_nodes[5][4] =
{
	{0, 1, 4, 99}, // Face 0
	{1, 2, 4, 99}, // Face 1
	{2, 3, 4, 99}, // Face 2
	{3, 0, 4, 99}, // Face 3
	{0, 3, 2,  1}  // Face 4
};

const unsigned PyramidRule5::edge_nodes[8][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}, // Edge 3
	{0, 4}, // Edge 4
	{1, 4}, // Edge 5
	{2, 4}, // Edge 6
	{3, 4}  // Edge 7
};

const unsigned PyramidRule5::n_face_nodes[5] = { 3, 3, 3, 3, 4 };

unsigned PyramidRule5::getNFaceNodes(unsigned i)
{
	if (i<5)
		return n_face_nodes[i];
	ERR("Error in MeshLib::Element::getNFaceNodes() - Index %d does not exist.", i);
	return 0;
}

const Element* PyramidRule5::getFace(const Element* e, unsigned i)
{
	if (i<e->getNFaces())
	{
		unsigned nFaceNodes (e->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));

		if (i<4)
			return new Tri(nodes);
		else
			return new Quad(nodes);
	}
	ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
	return nullptr;
}

double PyramidRule5::computeVolume(Node const* const* _nodes)
{
	return GeoLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[4]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[2]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords(), _nodes[4]->getCoords());
}

bool PyramidRule5::isPntInElement(Node const* const* _nodes, MathLib::Point3d const& pnt, double eps)
{
	return (GeoLib::isPointInTetrahedron(pnt, *_nodes[0], *_nodes[1], *_nodes[2], *_nodes[4], eps) ||
			GeoLib::isPointInTetrahedron(pnt, *_nodes[0], *_nodes[2], *_nodes[3], *_nodes[4], eps));
}

unsigned PyramidRule5::identifyFace(Node const* const* _nodes, Node* nodes[3])
{
	for (unsigned i=0; i<5; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<4; j++)
			for (unsigned k=0; k<3; k++)
				if (face_nodes[i][j] != 99 && _nodes[face_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

ElementErrorCode PyramidRule5::validate(const Element* e)
{
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();

	const MeshLib::Quad* base (dynamic_cast<const MeshLib::Quad*>(e->getFace(4)));
	if (base)
	{
		error_code |= base->validate();
		error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
	}
	else
		error_code.set(ElementErrorFlag::NodeOrder);
	delete base;

	return error_code;
}

} // end namespace MeshLib
