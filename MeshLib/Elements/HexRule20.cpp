/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HexRule20.h"

#include "logog/include/logog.hpp"

#include "Quad.h"

namespace MeshLib
{

const unsigned HexRule20::n_all_nodes;

const CellType HexRule20::cell_type;

//TODO
const unsigned HexRule20::_face_nodes[6][8] =
{
	{0, 3, 2, 1, 0, 0, 0, 0}, // Face 0
	{0, 1, 5, 4, 0, 0, 0, 0}, // Face 1
	{1, 2, 6, 5, 0, 0, 0, 0}, // Face 2
	{2, 3, 7, 6, 0, 0, 0, 0}, // Face 3
	{3, 0, 4, 7, 0, 0, 0, 0}, // Face 4
	{4, 5, 6, 7, 0, 0, 0, 0}  // Face 5
};

//TODO
const unsigned HexRule20::_edge_nodes[12][3] =
{
	{0, 1, 0}, // Edge 0
	{1, 2, 0}, // Edge 1
	{2, 3, 0}, // Edge 2
	{0, 3, 0}, // Edge 3
	{4, 5, 0}, // Edge 4
	{5, 6, 0}, // Edge 5
	{6, 7, 0}, // Edge 6
	{4, 7, 0}, // Edge 7
	{0, 4, 0}, // Edge 8
	{1, 5, 0}, // Edge 9
	{2, 6, 0}, // Edge 10
	{3, 7, 0}  // Edge 11
};

const Element* HexRule20::getFace(Node const* const* _nodes, unsigned i)
{
	if (i < n_faces)
	{
		unsigned nFaceNodes (getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = const_cast<Node*>(_nodes[_face_nodes[i][j]]);
		return new Quad8(nodes);
	}
	ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
	return NULL;
}


} // end namespace MeshLib
