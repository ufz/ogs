/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Pyramid.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef PYRAMID_H_
#define PYRAMID_H_

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
class Pyramid : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Pyramid(Node* nodes[5], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Pyramid(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, unsigned value = 0);

	/// Copy constructor
	Pyramid(const Pyramid &pyramid);

	/// Destructor
	virtual ~Pyramid();

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
	virtual unsigned getNNodes() const { return 5; };

	/**
	 * Method returns the type of the element. In this case PYRAMID will be returned.
	 * @return MshElemType::PYRAMID
	 */
	virtual MshElemType::type getType() const { return MshElemType::PYRAMID; }

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the
	 * Pyramid instance employing the copy constructor of class Pyramid.
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
	inline Node const* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

	static const unsigned _face_nodes[5][4];
	static const unsigned _edge_nodes[8][2];
	static const unsigned _n_face_nodes[5];

}; /* class */

} /* namespace */

#endif /* PYRAMID_H_ */

