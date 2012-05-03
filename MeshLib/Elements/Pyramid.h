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
 * @code
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
 *
 * @endcode
 */
class Pyramid : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Pyramid(Node* nodes[5], size_t value = 0);

	/// Constructor using single mesh nodes.
	Pyramid(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, size_t value = 0);

	/// Copy constructor
	Pyramid(const Pyramid &pyramid);

	/// Destructor
	virtual ~Pyramid();

	/// Get the number of nodes for this element.
	size_t getNNodes() const { return 5; };

protected:
	/// Calculates the volume of a prism by subdividing it into two tetrahedra.
	double calcVolume();

}; /* class */

} /* namespace */

#endif /* PYRAMID_H_ */

