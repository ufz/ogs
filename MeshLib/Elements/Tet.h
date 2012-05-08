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
 * @code
 *
 *  Tet:  3
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
class Tet : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Tet(Node* nodes[4], size_t value = 0);

	/// Constructor using single mesh nodes.
	Tet(Node* n0, Node* n1, Node* n2, Node* n3, size_t value = 0);

	/// Copy constructor
	Tet(const Tet &tet);

	/// Destructor
	virtual ~Tet();

	/// Get the number of edges for this element.
	size_t getNEdges() const { return 6; };
	
	/// Get the number of faces for this element.
	size_t getNFaces() const { return 4; };

	/// Get the number of neighbors for this element.
	size_t getNNeighbors() const { return 4; };

	/// Get the number of nodes for this element.
	virtual size_t getNNodes() const { return 4; };

protected:
	/// Constructor without nodes (for use of derived classes)
	Tet(size_t value = 0);

	/// Calculates the volume of a tetrahedron via the determinant of the matrix given by its four points.
	double calcVolume();

}; /* class */

} /* namespace */

#endif /* TET_H_ */

