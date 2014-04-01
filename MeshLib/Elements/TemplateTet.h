/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the TemplateTet class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATETET_H_
#define TEMPLATETET_H_

#include <array>
#include "MeshEnums.h"
#include "Cell.h"

namespace MeshLib {

/**
 * This class represents a 3d tetrahedron element. The following sketch shows the node and edge numbering.
 * @anchor TetrahedronNodeAndEdgeNumbering
 * @code
 *          3
 *         /|\
 *        / | \
 *      3/  |  \5
 *      /   |4  \
 *     /    |    \
 *    0.....|.....2
 *     \    |  2 /
 *      \   |   /
 *      0\  |  /1
 *        \ | /
 *         \|/
 *          1
 *
 * @endcode
 */
template <unsigned NNODES, CellType CELLTETTYPE>
class TemplateTet : public Cell
{
public:
	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes;

	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes;

	/// Constructor with an array of mesh nodes.
	TemplateTet(Node* nodes[NNODES], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs a tetrahedron from array of Node pointers.
	TemplateTet(std::array<Node*, NNODES> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Copy constructor
	TemplateTet(const TemplateTet &tet);

	/// Destructor
	virtual ~TemplateTet();

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const;

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 6; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { (void)i; return 3; };

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 4; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 4; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(bool all = false) const
	{
		return all ? n_all_nodes : n_base_nodes;
	}

	/**
	 * Method returns the type of the element. In this case TETRAHEDRON will be returned.
	 * @return MeshElemType::TETRAHEDRON
	 */
	virtual MeshElemType getGeomType() const { return MeshElemType::TETRAHEDRON; }

	/**
	 * Get the type of the element in context of the finite element method.
	 * @return a value of the enum CellType
	 */
	virtual CellType getCellType() const { return CELLTETTYPE; }

	/// Returns true if these two indeces form an edge and false otherwise
	bool isEdge(unsigned i, unsigned j) const;

	/**
	 * Tests if the element is geometrically valid.
	 * @param check_zero_volume indicates if volume == 0 should be checked
	 */
	virtual ElementErrorCode validate() const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateTet instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

protected:
	/// Calculates the volume of a tetrahedron via the determinant of the matrix given by its four points.
	double computeVolume();

	/**
	 * Return a specific edge node, see @ref TetrahedronNodeAndEdgeNumbering for numbering
	 * @param edge_id the id/number of the edge, have to be an integer value in the interval [0,6]
	 * @param node_id the id of the node within the edge (either 0 or 1)
	 * @return a pointer to the internal Node
	 */
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const {
		return _nodes[_edge_nodes[edge_id][node_id]];
	}

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;


	static const unsigned _face_nodes[4][3];
	static const unsigned _edge_nodes[6][2];

}; /* class */

} /* namespace */

#include "TemplateTet-impl.h"

#endif /* TEMPLATETET_H_ */

