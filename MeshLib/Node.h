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
	Node(double const*const coords, size_t id = std::numeric_limits<size_t>::max());
	Node(double x, double y, double z, size_t id = std::numeric_limits<size_t>::max());
	Node(const Node &node);
	~Node();

private:
	std::vector<Element*> _elements;

}; /* class */

} /* namespace */

#endif /* NODE_H_ */
