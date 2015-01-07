/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef HEXRULE8_H_
#define HEXRULE8_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "CellRule.h"

namespace MeshLib
{

/**
 * A 8-nodes Hexahedron Element.
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
class HexRule8 : public CellRule
{
public:
	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes = 8u;

	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = 8u;

	/// Constant: The geometric type of the element
	static const MeshElemType mesh_elem_type = MeshElemType::HEXAHEDRON;

	/// Constant: The FEM type of the element
	static const CellType cell_type = CellType::HEX8;

	/// Constant: The number of faces
	static const unsigned n_faces = 6;

	/// Constant: The number of edges
	static const unsigned n_edges = 12;

	/// Constant: The number of neighbors
	static const unsigned n_neighbors = 6;

	/// Constant: Local node index table for faces
	static const unsigned _face_nodes[6][4];

	/// Constant: Local node index table for edge
	static const unsigned _edge_nodes[12][2];

	/// Get the number of nodes for face i.
	static unsigned getNFaceNodes(unsigned /*i*/) { return 4; }

	/// Returns the i-th face of the element.
	static const Element* getFace(Node const* const* nodes, unsigned i);

	/**
	 * Checks if a point is inside the element.
	 * @param pnt a 3D GeoLib::Point object
	 * @param eps tolerance for numerical algorithm used or computing the property
	 * @return true if the point is not outside the element, false otherwise
	 */
	static bool isPntInElement(Node const* const* _nodes, GeoLib::Point const& pnt, double eps = std::numeric_limits<double>::epsilon());

	/**
	 * Tests if the element is geometrically valid.
	 * @param check_zero_volume indicates if volume == 0 should be checked
	 */
	static ElementErrorCode validate(const Element* e);

	/// Returns the ID of a face given an array of nodes.
	static unsigned identifyFace(Node const* const*, Node* nodes[3]);

	/// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
	static double computeVolume(Node const* const* _nodes);

}; /* class */

} /* namespace */

#endif /* HEXRULE_H_ */

