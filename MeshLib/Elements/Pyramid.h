/**
 * Pyramid.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef PYRAMID_H_
#define PYRAMID_H_

#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Pyramid Element.
 *
 *  Pyramid:   o 4
 *           //|\
 *          // | \
 *         //  |  \
 *      3 o/...|...o 2
 *       ./    |  /
 *      ./     | /
 *     ./      |/
 *    o--------o
 *    0        1
 */
class Pyramid : public Cell
{
public:
	Pyramid(Node* nodes[5], size_t value = 0);
	Pyramid(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, size_t value = 0);
	Pyramid(const Pyramid &pyramid);
	virtual ~Pyramid();

	size_t getNNodes() const { return 5; };

protected:
	/// Calculates the volume of a prism by subdividing it into two tetrahedra.
	double calcVolume();

}; /* class */

} /* namespace */

#endif /* PYRAMID_H_ */
