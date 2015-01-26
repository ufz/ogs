/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TRIRULE3_H_
#define TRIRULE3_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "FaceRule.h"
#include "EdgeReturn.h"

namespace MeshLib
{

/**
 * This class represents a 2d triangle element. The following sketch shows the node and edge numbering.
 * @anchor TriNodeAndEdgeNumbering
 * @code
 *
 *          2
 *          o
 *         / \
 *        /   \
 *      2/     \1
 *      /       \
 *     /         \
 *    0-----------1
 *          0
 *
 * @endcode
 */
class TriRule3 : public FaceRule
{
public:
	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes = 3u;

	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = 3u;

	/// Constant: The geometric type of the element
	static const MeshElemType mesh_elem_type = MeshElemType::TRIANGLE;

	/// Constant: The FEM type of the element
	static const CellType cell_type = CellType::TRI3;

	/// Constant: The number of faces
	static const unsigned n_faces = 0;

	/// Constant: The number of edges
	static const unsigned n_edges = 3;

	/// Constant: The number of neighbors
	static const unsigned n_neighbors = 3;

	/// Constant: Local node index table for edge
	static const unsigned edge_nodes[3][2];

	/// Returns the i-th edge of the element.
	typedef LinearEdgeReturn EdgeReturn;

	/**
	 * Checks if a point is inside the element.
	 * @param pnt a 3D GeoLib::Point object
	 * @param eps tolerance for numerical algorithm used or computing the property
	 * @return true if the point is not outside the element, false otherwise
	 */
	static bool isPntInElement(Node const* const* _nodes, MathLib::Point3d const& pnt, double eps);

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

#endif /* TRIRULE3_H_ */

