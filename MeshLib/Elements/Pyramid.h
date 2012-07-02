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
 * A 3d Pyramid Element.
 * @code
 *
 *  Pyramid:   o 4
 *           //|\
 *          // | \
 *         //  |  \
 *      3 o/...|...o 2
 *       ./    |  /
 *      ./     | /
 *     ./      |/
 *    o--------o
 *    0        1
 *
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

protected:
	/// Calculates the volume of a prism by subdividing it into two tetrahedra.
	double computeVolume();

	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	static const unsigned _face_nodes[5][4];
	static const unsigned _edge_nodes[8][2];
	static const unsigned _n_face_nodes[5];

}; /* class */

} /* namespace */

#endif /* PYRAMID_H_ */

