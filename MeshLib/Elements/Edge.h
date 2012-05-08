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
	Edge(Node* nodes[2], size_t value = 0);

	/// Constructor using single nodes
	Edge(Node* n1, Node* n2, size_t value = 0);

	/// Copy constructor
	Edge(const Edge &edge);

	/// Destructor
	virtual ~Edge();

	/// 1D elements have no edges
	size_t getNEdges() const { return 0; };

	/// Get the number of faces for this element.
	size_t getNFaces() const { return 0; };

	/// Get the length of this 1d element.
	double getLength() const { return _length; };

	/// Get dimension of the mesh element.
	size_t getDimension() const { return 1; };

	/// Get the number of neighbors for this element.
	size_t getNNeighbors() const { return 0; };

	/// Get the number of nodes for this element.
	virtual size_t getNNodes() const { return 2; };


protected:
	/// Calculate the length of this 1d element.
	double calcLength();

	double _length;

}; /* class */

} /* namespace */

#endif /* EDGE_H_ */

