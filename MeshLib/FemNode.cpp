/**
 * FemNode.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "FemNode.h"

namespace MeshLib {

FemNode::FemNode(double const*const coords, size_t id)
	: Node(coords, id)
{
}

FemNode::FemNode(double x, double y, double z, size_t id)
	: Node(x, y, z, id)
{
}

FemNode::FemNode(const Node &node)
	: Node(node.getData(), node.getID())
{
}

FemNode::FemNode(const FemNode &node)
	: Node(node.getData(), node.getID())
{
}

FemNode::~FemNode()
{
}

}