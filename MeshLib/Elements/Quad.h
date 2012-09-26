/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
 * This class represents a 2d quadliteral element. The following sketch shows the node and edge numbering.
 * @anchor QuadNodeAndEdgeNumbering
 * @code
 *              2
 *        3-----------2
 *        |           |
 *        |           |
 *       3|           |1
 *        |           |
 *        |           |
 *        0-----------1
 *              0
 * @endcode
 */
template <unsigned ORDER, unsigned NNODES>
class TemplateQuad : public Face
{
public:
	/// Constructor with an array of mesh nodes.
	TemplateQuad(Node* nodes[NNODES], unsigned value = 0);

	/// Copy constructor
	TemplateQuad(const TemplateQuad &quad);

	/// Destructor
	virtual ~TemplateQuad();

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 4; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 4; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(unsigned order = 1) const
	{
		return order == ORDER ? NNODES : 4;
	}

	/**
	 * Method returns the type of the element. In this case QUAD will be returned.
	 * @return MshElemType::QUAD
	 */
	virtual MshElemType::type getType() const { return MshElemType::QUAD; }

	/// Returns true if these two indeces form an edge and false otherwise
	bool isEdge(unsigned i, unsigned j) const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateQuad instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

	/**
	 * This method should be called after at least two nodes of the quad
	 * element are collapsed. As a consequence of the node collapsing an edge
	 * of the quad will be collapsed. If one of the edges (see
	 * sketch @ref PyramidNodeAndEdgeNumbering) is collapsed we obtain a
	 * triangle. In this case the method will create the appropriate
	 * object of class Tri.
	 * @return a Tri object or NULL
	 */
	virtual Element* reviseElement() const;

protected:
	/// Calculates the area of a convex quadliteral by dividing it into two triangles.
	double computeVolume();

protected:
	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

	static const unsigned _edge_nodes[4][2];

}; /* class */

template <unsigned ORDER, unsigned NNODES>
const unsigned TemplateQuad<ORDER,NNODES>::_edge_nodes[4][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}  // Edge 3
};

typedef TemplateQuad<1,4> Quad;

} /* namespace */

#include "Quad.hpp"

#endif /* QUAD_H_ */

