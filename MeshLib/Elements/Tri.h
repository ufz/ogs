/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file Tri.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef TRI_H_
#define TRI_H_

#include "Face.h"

namespace MeshLib {

/**
 * A 2d Triangle Element.
 * @code
 *
 *  Tri:    2
 *          o
 *         / \
 *        /   \
 *       /     \
 *      /       \
 *     /         \
 *    o-----------o
 *    0           1
 *
 * @endcode
 */
class Tri : public Face
{
public:
	/// Constructor with an array of mesh nodes.
	Tri(Node* nodes[3], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Tri(Node* n0, Node* n1, Node* n2, unsigned value = 0);

	/// Copy constructor
	Tri(const Tri &tri);

	/// Destructor
	virtual ~Tri();

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 3; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 3; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const { return 3; };

protected:
	/// Calculates the area of the triangle by returning half of the area of the corresponding parallelogram.
	double computeArea();

	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	static const unsigned _edge_nodes[3][2];

}; /* class */

} /* namespace */

#endif /* TRI_H_ */

