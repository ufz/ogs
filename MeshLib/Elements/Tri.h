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
 *
 *   Tri:   2
 *          o
 *         / \
 *        /   \
 *       /     \
 *      /       \
 *     /         \
 *    o-----------o
 *    0           1
 *
 */
class Tri : public Face
{
public:
	Tri(Node* nodes[3], size_t value = 0);
	Tri(Node* n0, Node* n1, Node* n2, size_t value = 0);
	Tri(const Tri &tri);
	virtual ~Tri();

	size_t getNNodes() const { return 3; };

protected:
	/// Calculates the area of the triangle by returning half of the area of the corresponding parallelogram.
	double calcArea();

}; /* class */

} /* namespace */

#endif /* TRI_H_ */
