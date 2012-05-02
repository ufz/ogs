/**
 * Node.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Node.h"

namespace MeshLib {

Node::Node(const double coords[3], size_t id)
	: GEOLIB::PointWithID(coords, id)
{
}

Node::Node(double x, double y, double z, size_t id)
	: GEOLIB::PointWithID(x, y, z, id)
{
}

Node::Node(const Node &node)
	: GEOLIB::PointWithID(node.getData(), node.getID())
{
}

Node::~Node()
{
}

}