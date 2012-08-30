/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Quad.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef QUAD_H_
#define QUAD_H_

#include "Face.h"

namespace MeshLib {

/**
 * A 2d Quadliteral Element.
 * @code
 *
 *        3           2
 *  Quad: o-----------o
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *        o-----------o
 *        0           1
 *
 * @endcode
 */
class Quad : public Face
{
public:
	/// Constructor with an array of mesh nodes.
	Quad(Node* nodes[4], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Quad(Node* n0, Node* n1, Node* n2, Node* n3, unsigned value = 0);

	/// Copy constructor
	Quad(const Quad &quad);

	/// Destructor
	virtual ~Quad();

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 4; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 4; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const { return 4; };

	/**
	 * Method returns the type of the element. In this case QUAD will be returned.
	 * @return MshElemType::QUAD
	 */
	virtual MshElemType::type getType() const { return MshElemType::QUAD; }

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the Quad instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

	/// Calculates the area of a convex quadliteral by dividing it into two triangles.
	double computeVolume();

protected:
	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

	static const unsigned _edge_nodes[4][2];

}; /* class */

} /* namespace */

#endif /* QUAD_H_ */

