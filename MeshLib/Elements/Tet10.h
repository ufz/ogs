/**
 * Tet10.h
 *
 *      Date: 2012/05/03
 *      Author: KR
 */

#ifndef TET10_H_
#define TET10_H_

#include "Tet.h"
#include "FemElem.h"

namespace MeshLib {

/**
 * A 3d Tetrahedron Finite Element with 10 Nodes.
 * @code
 *
 *  Tet10:    3
 *            o
 *           /|\
 *          / | \
 *      7  /  |  \9
 *        o   |   o
 *       /    |8   \
 *      /     o     \
 *     /    6 |      \
 *  0 o.....o.|.......o 2
 *     \      |      /
 *      \     |     /
 *       \    |    /
 *      4 o   |   o 5
 *         \  |  /
 *          \ | /
 *           \|/
 *            o
 *            1
 * @endcode
 */
class Tet10 : public Tet, public FemElem
{
public:
	/// Constructor with an array of mesh nodes.
	Tet10(Node* nodes[10], size_t value = 0);

	/// Constructor using a simple Tetrahedron
	Tet10(const Tet &tet);

	/// Copy constructor
	Tet10(const Tet10 &tet);

	/// Destructor
	virtual ~Tet10();

	/// Get the number of nodes for this element.
	size_t getNNodes() const { return 10; };

protected:
	/// Calculates the volume of a tetrahedron via the determinant of the matrix given by its four points.
	void calcCoG();

}; /* class */

} /* namespace */

#endif /* TET10_H_ */

