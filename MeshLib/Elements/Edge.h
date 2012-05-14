/**
 * Edge.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "Element.h"

namespace MeshLib {

class Node;

/**
 * A 1d Edge or Line Element.
 * @code
 *
 *  Edge: o--------o
 *        0        1
 *
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


protected:
	/// Calculate the length of this 1d element.
	double computeLength();

	/// 1D elements have no edges.
	Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { (void)edge_id; (void)node_id; return NULL; };

	/// 1D elements have no faces.
	Node* getFaceNode(unsigned face_id, unsigned node_id) const { (void)face_id; (void)node_id; return NULL; };

	double _length;

}; /* class */

} /* namespace */

#endif /* EDGE_H_ */

