/**
 * \file FemNode.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "FemNode.h"

namespace MeshLib {

FemNode::FemNode(double const*const coords, unsigned id)
	: Node(coords, id)
{
}

FemNode::FemNode(double x, double y, double z, unsigned id)
	: Node(x, y, z, id)
{
}

FemNode::FemNode(const Node &node)
	: Node(node.getCoords(), node.getID())
{
}

FemNode::FemNode(const FemNode &node)
	: Node(node.getCoords(), node.getID())
{
}

FemNode::~FemNode()
{
}

}

