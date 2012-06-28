/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file Hex.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef HEX_H_
#define HEX_H_

#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Hexahedron Element.
 * @code
 *
 *  Hex:  7        6
 *        o--------o
 *       /:       /|
 *      / :      / |
 *   4 /  :   5 /  |
 *    o--------o   |
 *    |   o....|...o 2
 *    |  .3    |  /
 *    | .      | /
 *    |.       |/
 *    o--------o
 *    0        1
 *
 * @endcode
 */
class Hex : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Hex(Node* nodes[8], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Hex(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, Node* n6, Node* n7, unsigned value);

	/// Copy constructor
	Hex(const Hex &hex);

	/// Destructor
	virtual ~Hex();

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
	virtual unsigned getNNodes() const { return 8; };

protected:
	/// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
	double computeVolume();

	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	static const unsigned _face_nodes[6][4];
	static const unsigned _edge_nodes[12][2];

}; /* class */

} /* namespace */

#endif /* HEX_H_ */

