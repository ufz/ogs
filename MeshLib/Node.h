/**
 * Node.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef NODE_H_
#define NODE_H_

#include <cstdlib>
#include <limits>
#include <vector>
#include "PointWithID.h"

namespace MeshLib {

class Element;

/**
 * A mesh node with coordinates in 3D space.
 */
class Node : public GEOLIB::PointWithID
{
public:
	/// Constructor using a coordinate array
	Node(const double coords[3], size_t id = std::numeric_limits<size_t>::max());

	/// Constructor using single coordinates
	Node(double x, double y, double z, size_t id = std::numeric_limits<size_t>::max());

	/// Copy constructor
	Node(const Node &node);

	/// Destructor
	~Node();

private:
	std::vector<Element*> _elements;

}; /* class */

} /* namespace */

#endif /* NODE_H_ */
