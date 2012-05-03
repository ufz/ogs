/**
 * Prism.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef PRISM_H_
#define PRISM_H_

#include "Cell.h"

namespace MeshLib {

/**
 * A 3d Prism Element.
 * @code
 *
 *  Prism:   5
 *           o
 *          /:\
 *         / : \
 *        /  o  \
 *     3 o-------o 4
 *       | . 2 . |
 *       |.     .|
 *       o-------o
 *       0       1
 *
 * @endcode
 */
class Prism : public Cell
{
public:
	/// Constructor with an array of mesh nodes.
	Prism(Node* nodes[6], size_t value = 0);

	/// Constructor using single mesh nodes.
	Prism(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, size_t value = 0);

	/// Copy constructor
	Prism(const Prism &prism);

	/// Destructor
	virtual ~Prism();

	/// Get the number of nodes for this element.
	size_t getNNodes() const { return 6; };

protected:
	/// Calculates the volume of a prism by subdividing it into three tetrahedra.
	double calcVolume();

}; /* class */

} /* namespace */

#endif /* PRISM_H_ */

