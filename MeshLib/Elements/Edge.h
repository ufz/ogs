/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Edge.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <limits>

#include "Element.h"

namespace MeshLib {

class Node;

/**
 * A 1d Edge or Line Element.
 * @code
 *  0--------1
 * @endcode
 */
class Edge : public Element
{
public:
	/// Constructor with an array of mesh nodes.
	Edge(Node* nodes[2], unsigned value = 0);

	/// Constructor using single nodes
	Edge(Node* n1, Node* n2, unsigned value = 0);

	/// Copy constructor
	Edge(const Edge &edge);

	/// Destructor
	virtual ~Edge();

	/// Returns the length, area or volume of a 1D, 2D or 3D element
	double getContent() const { return _length; };

	/// Returns the edge i of the element.
	const Element* getEdge(unsigned i) const { (void)i; return NULL; };

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const { (void)i; return NULL; };

	/// Compute the minimum and maximum squared edge length for this element
	void computeSqrEdgeLengthRange(double &min, double &max) const { min = _length; max = _length; };

	/// 1D elements have no edges
	unsigned getNEdges() const { return 0; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { (void)i; return 0; };

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 0; };

	/// Get the length of this 1d element.
	double getLength() const { return _length; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return 1; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 0; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const { return 2; };

	virtual MshElemType::type getType() const { return MshElemType::EDGE; }

	virtual Element* clone() const;

	/**
	 * If for instance two nodes of the element are collapsed the Edge element disappears.
	 * @return NULL
	 */
	virtual Element* reviseElement() const;

protected:
	double computeVolume();
	/// 1D elements have no edges.
	Node const* getEdgeNode(unsigned edge_id, unsigned node_id) const { (void)edge_id; (void)node_id; return NULL; };

	/// 1D elements have no faces.
	Node* getFaceNode(unsigned face_id, unsigned node_id) const { (void)face_id; (void)node_id; return NULL; };

	/// Returns the ID of a face given an array of nodes (but is not applicable for edges!).
	unsigned identifyFace(Node* [3]/*nodes[3]*/) const { return std::numeric_limits<unsigned>::max(); };

	double _length;

}; /* class */

} /* namespace */

#endif /* EDGE_H_ */

