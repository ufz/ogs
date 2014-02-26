/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the TemplatePyramid class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEPYRAMID_H_
#define TEMPLATEPYRAMID_H_

#include <array>
#include "MeshEnums.h"
#include "Cell.h"

namespace MeshLib {

/**
 * This class represents a 3d pyramid element. The following sketch shows the node and edge numbering.
 * @anchor PyramidNodeAndEdgeNumbering
 * @code
 *
 *               4
 *             //|\
 *            // | \
 *          7//  |  \6
 *          //   |5  \
 *         //    |    \
 *        3/.... |.....2
 *       ./      |  2 /
 *      ./4      |   /
 *    3./        |  /1
 *    ./         | /
 *   ./          |/
 *  0------------1
 *        0
 * @endcode
 */
template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
class TemplatePyramid : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	TemplatePyramid(Node* nodes[NNODES], unsigned value = 0);

	/// Constructs a pyramid from array of Node pointers.
	TemplatePyramid(std::array<Node*, NNODES> const& nodes, unsigned value = 0);

	/// Copy constructor
	TemplatePyramid(const TemplatePyramid &pyramid);

	/// Destructor
	virtual ~TemplatePyramid();

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const;

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 8; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const;

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 5; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 5; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(bool all = false) const
	{
		return all ? NNODES : 5;
	}

	/**
	 * Method returns the type of the element. In this case PYRAMID will be returned.
	 * @return MeshElemType::PYRAMID
	 */
	virtual MeshElemType getGeomType() const { return MeshElemType::PYRAMID; }

	/**
	 * Get the type of the element in context of the finite element method.
	 * @return a value of the enum CellType
	 */
	virtual CellType getCellType() const { return CELLPYRAMIDTYPE; }

	/// Returns true if these two indeces form an edge and false otherwise
	bool isEdge(unsigned i, unsigned j) const;

	/**
	 * Tests if the element is geometrically valid.
	 * @param check_zero_volume indicates if volume == 0 should be checked
	 */
	virtual ElementErrorCode validate() const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the
	 * TemplatePyramid instance employing the copy constructor of class TemplatePyramid.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

	/**
	 * This method should be called after at least two nodes of the pyramid
	 * element are collapsed. As a consequence of the node collapsing an edge
	 * of the pyramid will be collapsed. If one of the edges 0, 1, 2 or 3 (see
	 * sketch @ref PyramidNodeAndEdgeNumbering) is collapsed we obtain a
	 * tetrahedron. In this case the method will create the appropriate
	 * object of class Tetrahedron.
	 * @return a Tetrahedron object or NULL
	 */
	virtual Element* reviseElement() const;

protected:
	/// Calculates the volume of a prism by subdividing it into two tetrahedra.
	double computeVolume();

	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

	static const unsigned _face_nodes[5][4];
	static const unsigned _edge_nodes[8][2];
	static const unsigned _n_face_nodes[5];

}; /* class */

} /* namespace */

#include "TemplatePyramid.tpp"

#endif /* TEMPLATEPYRAMID_H_ */

