/**
 * Hex.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef HEX_H_
#define HEX_H_

#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Hexahedron Element.
 * @code
 *
 *  Hex:  7        6
 *        o--------o
 *       /:       /|
 *      / :      / |
 *   4 /  :   5 /  |
 *    o--------o   |
 *    |   o....|...o 2
 *    |  .3    |  /
 *    | .      | /
 *    |.       |/
 *    o--------o
 *    0        1
 *
 * @endcode
 */
class Hex : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Hex(Node* nodes[8], unsigned value = 0);

	/// Constructor using single mesh nodes.
	Hex(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, Node* n6, Node* n7, unsigned value);

	/// Copy constructor
	Hex(const Hex &hex);

	/// Destructor
	virtual ~Hex();

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 12; };

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 6; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 6; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const { return 8; };

protected:
	/// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
	double calcVolume();

}; /* class */

} /* namespace */

#endif /* HEX_H_ */

