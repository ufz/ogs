/**
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEHEX20_H_
#define TEMPLATEHEX20_H_

#include <array>
#include "MeshEnums.h"
#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Hexahedron Element.
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
template <unsigned NNODES, CellType CELLHEXTYPE>
class TemplateHex20 : public Cell
{
public:
	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = NNODES;

	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes = 8u;

	/// Constructor with an array of mesh nodes.
	TemplateHex20(Node* nodes[NNODES], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs a hex from array of Node pointers.
	TemplateHex20(std::array<Node*, NNODES> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Copy constructor
	TemplateHex20(const TemplateHex20 &hex);

	/// Destructor
	virtual ~TemplateHex20();

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const;

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 12; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { (void)i; return 4; };

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 6; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 6; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const
	{
		return n_all_nodes;
	}

	/// Get the number of linear nodes for this element.
	virtual unsigned getNLinearNodes() const
	{
		return n_base_nodes;
	}

	/**
	 * Method returns the type of the element. In this case HEXAHEDRON will be returned.
	 * @return MeshElemType::HEXAHEDRON
	 */
	virtual MeshElemType getGeomType() const { return MeshElemType::HEXAHEDRON; }

	/**
	 * Method returns the FEM type of the element.
	 * @return
	 */
	virtual CellType getCellType() const { return CELLHEXTYPE; }

	/// Returns true if these two indices form an edge and false otherwise
	bool isEdge(unsigned i, unsigned j) const;

	/**
	 * Tests if the element is geometrically valid.
	 * @param check_zero_volume indicates if volume == 0 should be checked
	 */
	virtual ElementErrorCode validate() const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the Hex instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

protected:
	/// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
	double computeVolume();

	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	static const unsigned _face_nodes[6][8];
	static const unsigned _edge_nodes[12][3];

}; /* class */

} /* namespace */

#include "TemplateHex20-impl.h"

#endif /* TEMPLATEHEX_H_ */

