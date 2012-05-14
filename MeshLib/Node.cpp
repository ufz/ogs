/**
 * Node.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Node.h"

namespace MeshLib {

Node::Node(const double coords[3], unsigned id)
	: GeoLib::PointWithID(coords, id)
{
}

Node::Node(double x, double y, double z, unsigned id)
	: GeoLib::PointWithID(x, y, z, id)
{
}

Node::Node(const Node &node)
	: GeoLib::PointWithID(node.getCoords(), node.getID())
{
}

Node::~Node()
{
}

}

