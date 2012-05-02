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
 *
 *  Edge: o--------o
 *        0        1
 */
class Edge : public Element
{
public:
	Edge(Node* nodes[2], size_t value = 0);
	Edge(Node* n1, Node* n2, size_t value = 0);
	Edge(const Edge &edge);
	virtual ~Edge();

	double getLength() const { return _length; };

	size_t getDimension() const { return 1; };

	size_t getNNodes() const { return 2; };


protected:
	double calcLength();

	double _length;

}; /* class */

} /* namespace */

#endif /* EDGE_H_ */
