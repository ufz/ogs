/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef HEXRULE20_H_
#define HEXRULE20_H_

#include "MeshLib/MeshEnums.h"
#include "HexRule8.h"

namespace MeshLib
{

/**
 * A 20-nodes Hexahedron Element.
 * @code
 *
 *  Hex:
 *                6
 *          7-----------6
 *         /:          /|
 *        / :         / |
 *      7/  :        /5 |
 *      / 11:       /   | 10
 *     /    : 4    /    |
 *    4-----------5     |
 *    |     :     | 2   |
 *    |     3.....|.....2
 *    |    .      |    /
 *  8 |   .       |9  /
 *    | 3.        |  / 1
 *    | .         | /
 *    |.          |/
 *    0-----------1
 *          0
 *
 * @endcode
 */
class HexRule20 : public HexRule8
{
public:
	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = 20u;

	/// Constant: The FEM type of the element
	static const CellType cell_type = CellType::HEX20;

	static const unsigned _face_nodes[6][8];
	static const unsigned _edge_nodes[12][3];

	/// Get the number of nodes for face i.
	static unsigned getNFaceNodes(unsigned /*i*/) { return 8; }

	/// Returns the i-th face of the element.
	static const Element* getFace(Node const* const* _nodes, unsigned i);

}; /* class */

} /* namespace */

#endif /* HEXRULE_H_ */

