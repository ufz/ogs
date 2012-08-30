/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Prism.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef PRISM_H_
#define PRISM_H_

#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Prism Element.
 * @code
 *
 *  Prism:   5
 *           o
 *          /:\
 *         / : \
 *        /  o  \
 *     3 o-------o 4
 *       | . 2 . |
 *       |.     .|
 *       o-------o
 *       0       1
 *
 * @endcode
 */
class Prism : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Prism(Node* nodes[6], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Prism(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, unsigned value = 0);

	/// Copy constructor
	Prism(const Prism &prism);

	/// Destructor
	virtual ~Prism();

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
	virtual unsigned getNNodes() const { return 6; };

	/**
	 * Method returns the type of the element. In this case PRISM will be returned.
	 * @return MshElemType::PRISM
	 */
	virtual MshElemType::type getType() const { return MshElemType::PRISM; }

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the
	 * Hex instance employing the copy constructor of class Prism.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

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

#endif /* PRISM_H_ */

