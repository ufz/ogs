/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tet.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef TET_H_
#define TET_H_

#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Tetrahedron Element.
 * @code
 *
 *  Tet:  3
 *        o
 *       /|\
 *      / | \
 *     /  |  \
 *  0 o...|...o 2
 *     \  |  /
 *      \ | /
 *       \|/
 *        o
 *        1
 *
 * @endcode
 */
class Tet : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Tet(Node* nodes[4], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Tet(Node* n0, Node* n1, Node* n2, Node* n3, unsigned value = 0);

	/// Copy constructor
	Tet(const Tet &tet);

	/// Destructor
	virtual ~Tet();

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const;

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 6; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { (void)i; return 3; };

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 4; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 4; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const { return 4; };

	/**
	 * Method returns the type of the element. In this case TETRAHEDRON will be returned.
	 * @return MshElemType::TETRAHEDRON
	 */
	virtual MshElemType::type getType() const { return MshElemType::TETRAHEDRON; }

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the Tet instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

protected:
	/// Constructor without nodes (for use of derived classes)
	Tet(unsigned value = 0);

	/// Calculates the volume of a tetrahedron via the determinant of the matrix given by its four points.
	double computeVolume();

	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	static const unsigned _face_nodes[4][3];
	static const unsigned _edge_nodes[6][2];

}; /* class */

} /* namespace */

#endif /* TET_H_ */

