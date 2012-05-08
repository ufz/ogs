/**
 * Tri.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef TRI_H_
#define TRI_H_

#include "Face.h"

namespace MeshLib {

/**
 * A 2d Triangle Element.
 * @code
 *
 *  Tri:    2
 *          o
 *         / \
 *        /   \
 *       /     \
 *      /       \
 *     /         \
 *    o-----------o
 *    0           1
 *
 * @endcode
 */
class Tri : public Face
{
public:
	/// Constructor with an array of mesh nodes.
	Tri(Node* nodes[3], size_t value = 0);

	/// Constructor using single mesh nodes.
	Tri(Node* n0, Node* n1, Node* n2, size_t value = 0);

	/// Copy constructor
	Tri(const Tri &tri);

	/// Destructor
	virtual ~Tri();

	/// Get the number of edges for this element.
	size_t getNEdges() const { return 3; };

	/// Get the number of neighbors for this element.
	size_t getNNeighbors() const { return 3; };

	/// Get the number of nodes for this element.
	virtual size_t getNNodes() const { return 3; };

protected:
	/// Calculates the area of the triangle by returning half of the area of the corresponding parallelogram.
	double calcArea();

}; /* class */

} /* namespace */

#endif /* TRI_H_ */

