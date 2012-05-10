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
 * @code
 *
 *        3           2
 *  Quad: o-----------o
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *        |           |
 *        o-----------o
 *        0           1
 *
 * @endcode
 */
class Quad : public Face
{
public:
	/// Constructor with an array of mesh nodes.
	Quad(Node* nodes[4], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Quad(Node* n0, Node* n1, Node* n2, Node* n3, unsigned value = 0);

	/// Copy constructor
	Quad(const Quad &quad);

	/// Destructor
	virtual ~Quad();

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 4; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 4; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const { return 4; };

protected:
	/// Calculates the area of a convex quadliteral by dividing it into two triangles.
	double calcArea();


}; /* class */

} /* namespace */

#endif /* QUAD_H_ */

