/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file TemplatePrism.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef TEMPLATEPRISM_H_
#define TEMPLATEPRISM_H_

#include "Cell.h"

namespace MeshLib {

/**
 * This class represents a 3d prism element. The following sketch shows the node and edge numbering.
 * @anchor PrismNodeAndEdgeNumbering
 * @code
 *            5
 *           / \
 *          / : \
 *        8/  :  \7
 *        /   :5  \
 *       /    :  6 \
 *      3-----------4
 *      |     :     |
 *      |     2     |
 *      |    . .    |
 *     3|   .   .   |4
 *      | 2.     .1 |
 *      | .       . |
 *      |.         .|
 *      0-----------1
 *            0
 *
 * @endcode
 */
template <unsigned NNODES>
class TemplatePrism : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	TemplatePrism(Node* nodes[6], unsigned value = 0);

	/// Copy constructor
	TemplatePrism(const TemplatePrism<NNODES> &prism);

	/// Destructor
	virtual ~TemplatePrism();

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const;

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 9; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const;

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 5; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 5; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(bool all = false) const
	{
		return all ? NNODES : 6;
	}

	/**
	 * Method returns the type of the element. In this case PRISM will be returned.
	 * @return MshElemType::PRISM
	 */
	virtual MshElemType::type getType() const { return MshElemType::PRISM; }

	/// Returns true if these two indeces form an edge and false otherwise
	bool isEdge(unsigned i, unsigned j) const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the
	 * Hex instance employing the copy constructor of class TemplatePrism.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

	/**
	 * This method should be called after at least two nodes of the prism
	 * element are collapsed. As a consequence of the node collapsing an edge
	 * of the prism will be collapsed. If one of the edges 3, 4 or 5 (see
	 * sketch @ref PrismNodeAndEdgeNumbering) is collapsed we obtain a
	 * pyramid. In this case the method will create the appropriate
	 * object of class Pyramid.
	 * @return a pyramid object or NULL
	 */
	virtual Element* reviseElement() const;

protected:
	/// Calculates the volume of a prism by subdividing it into three tetrahedra.
	double computeVolume();

	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

	static const unsigned _face_nodes[5][4];
	static const unsigned _edge_nodes[9][2];
	static const unsigned _n_face_nodes[5];

}; /* class */

} /* namespace */

#include "TemplatePrism.tpp"

#endif /* TEMPLATEPRISM_H_ */

