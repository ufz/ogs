/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HexRule8.h"

#include "logog/include/logog.hpp"

#include "GeoLib/AnalyticalGeometry.h"

#include "MeshLib/Node.h"
#include "Quad.h"

namespace MeshLib {

const unsigned HexRule8::n_all_nodes;

const unsigned HexRule8::n_base_nodes;

const unsigned HexRule8::_face_nodes[6][4] =
{
	{0, 3, 2, 1}, // Face 0
	{0, 1, 5, 4}, // Face 1
	{1, 2, 6, 5}, // Face 2
	{2, 3, 7, 6}, // Face 3
	{3, 0, 4, 7}, // Face 4
	{4, 5, 6, 7}  // Face 5
};

const unsigned HexRule8::_edge_nodes[12][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}, // Edge 3
	{4, 5}, // Edge 4
	{5, 6}, // Edge 5
	{6, 7}, // Edge 6
	{4, 7}, // Edge 7
	{0, 4}, // Edge 8
	{1, 5}, // Edge 9
	{2, 6}, // Edge 10
	{3, 7}  // Edge 11
};

const Element* HexRule8::getFace(Node const* const* _nodes, unsigned i)
{
	if (i < n_faces)
	{
		unsigned nFaceNodes (getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = const_cast<Node*>(_nodes[_face_nodes[i][j]]);
		return new Quad(nodes);
	}
	ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
	return NULL;
}

double HexRule8::computeVolume(Node const* const* _nodes)
{
	return GeoLib::calcTetrahedronVolume(_nodes[4]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[0]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[3]->getCoords(), _nodes[1]->getCoords(), _nodes[0]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[6]->getCoords(), _nodes[2]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[1]->getCoords(), _nodes[3]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[3]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords());
}

bool HexRule8::isPntInElement(Node const* const* _nodes, GeoLib::Point const& pnt, double eps)
{
	return (GeoLib::isPointInTetrahedron(pnt, *_nodes[4], *_nodes[7], *_nodes[5], *_nodes[0], eps) ||
			GeoLib::isPointInTetrahedron(pnt, *_nodes[5], *_nodes[3], *_nodes[1], *_nodes[0], eps) ||
			GeoLib::isPointInTetrahedron(pnt, *_nodes[5], *_nodes[7], *_nodes[3], *_nodes[0], eps) ||
			GeoLib::isPointInTetrahedron(pnt, *_nodes[5], *_nodes[7], *_nodes[6], *_nodes[2], eps) ||
			GeoLib::isPointInTetrahedron(pnt, *_nodes[1], *_nodes[3], *_nodes[5], *_nodes[2], eps) ||
			GeoLib::isPointInTetrahedron(pnt, *_nodes[3], *_nodes[7], *_nodes[5], *_nodes[2], eps));
}

unsigned HexRule8::identifyFace(Node const* const* _nodes, Node* nodes[3])
{
	for (unsigned i=0; i<6; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<4; j++)
			for (unsigned k=0; k<3; k++)
				if (_nodes[_face_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

ElementErrorCode HexRule8::validate(const Element* e)
{
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();
		
	for (unsigned i=0; i<6; ++i)
	{
		if (error_code.all())
			break;

		const MeshLib::Element* quad (e->getFace(i));
		error_code |= quad->validate();
		delete quad;
	}
	error_code[ElementErrorFlag::NodeOrder]  = !e->testElementNodeOrder();
	return error_code;
}

} // end namespace MeshLib
