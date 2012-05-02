/**
 * Quad.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef QUAD_H_
#define QUAD_H_

#include "Face.h"

namespace MeshLib {

/**
 * A 2d Quadliteral Element.
 *
 *        3           2
 * QUAD4: o-----------o
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *        o-----------o
 *        0           1
 */
class Quad : public Face
{
public:
	Quad(Node* nodes[4], size_t value = 0);
	Quad(Node* n0, Node* n1, Node* n2, Node* n3, size_t value = 0);
	Quad(const Quad &quad);
	virtual ~Quad();

	size_t getNNodes() const { return 4; };

protected:
	/// Calculates the area of a convex quadliteral by dividing it into two triangles.
	double calcArea();


}; /* class */

} /* namespace */

#endif /* QUAD_H_ */
