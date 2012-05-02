/**
 * Tet.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef TET_H_
#define TET_H_

#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Tetrahedron Element.
 * 
 * Tet:   3
 *        o
 *       /|\
 *      / | \
 *     /  |  \
 *  0 o...|...o 2
 *     \  |  /
 *      \ | /
 *       \|/
 *        o
 *        1
 */
class Tet : public Cell
{
public:
	Tet(Node* nodes[4], size_t value = 0);
	Tet(Node* n0, Node* n1, Node* n2, Node* n3, size_t value = 0);
	Tet(const Tet &tet);
	virtual ~Tet();

	size_t getNNodes() const { return 4; };

protected:
	/// Calculates the volume of a tetrahedron via the determinant of the matrix given by its four points.
	double calcVolume();

}; /* class */

} /* namespace */

#endif /* TET_H_ */
