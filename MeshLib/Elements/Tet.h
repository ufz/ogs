/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
 * This class represents a 3d tetrahedron element. The following sketch shows the node and edge numbering.
 * @anchor TetrahedronNodeAndEdgeNumbering
 * @code
 *          3
 *         /|\
 *        / | \
 *      3/  |  \5
 *      /   |4  \
 *     /    |    \
 *    0.....|.....2
 *     \    |  2 /
 *      \   |   /
 *      0\  |  /1
 *        \ | /
 *         \|/
 *          1
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

	/// Returns true if these two indeces form an edge and false otherwise
	bool isEdge(unsigned i, unsigned j) const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the Tet instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

	/**
	 * This method should be called after at least two nodes of the tetrahedron
	 * element are collapsed. As a consequence of the node collapsing an edge
	 * of the tetrahedron will be collapsed. If one edge is collapsed we obtain
	 * a triangle.
	 * @return a Triangle object or NULL
	 */
	virtual Element* reviseElement() const;

protected:
	/// Constructor without nodes (for use of derived classes)
	Tet(unsigned value = 0);

	/// Calculates the volume of a tetrahedron via the determinant of the matrix given by its four points.
	double computeVolume();

	/**
	 * Return a specific edge node, see @ref TetrahedronNodeAndEdgeNumbering for numbering
	 * @param edge_id the id/number of the edge, have to be an integer value in the interval [0,6]
	 * @param node_id the id of the node within the edge (either 0 or 1)
	 * @return a pointer to the internal Node
	 */
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;


	static const unsigned _face_nodes[4][3];
	static const unsigned _edge_nodes[6][2];

}; /* class */

} /* namespace */

#endif /* TET_H_ */

