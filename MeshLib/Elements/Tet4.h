/**
 * \file Tet4.h
 *
 * Created on 2012-05-03 by Karsten Rink
 */

#ifndef TET4_H_
#define TET4_H_

#include "Tet.h"
#include "FemElem.h"

namespace MeshLib {

/**
 * A 3d Tetrahedron Finite Element with 10 Nodes.
 * @code
 *
 *  Tet4: 3
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
 *
 * @endcode
 */
class Tet4 : public Tet, public FemElem
{
public:
	/// Constructor with an array of mesh nodes.
	Tet4(Node* nodes[4], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Tet4(Node* n0, Node* n1, Node* n2, Node* n3, unsigned value = 0);

	/// Constructor using a simple Tetrahedron
	Tet4(const Tet &tet);

	/// Copy constructor
	Tet4(const Tet4 &tet);

	/// Destructor
	virtual ~Tet4();

protected:
	/// Calculates the volume of a tetrahedron via the determinant of the matrix given by its four points.
	void calcCentroid();

}; /* class */

} /* namespace */

#endif /* TET4_H_ */

