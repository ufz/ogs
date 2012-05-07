/**
 * Edge.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Edge.h"
#include "Node.h"

#include "MathTools.h"

namespace MeshLib {

Edge::Edge(Node* nodes[2], size_t value)
	: Element(MshElemType::LINE, value)
{
	_nodes = nodes;
	this->_length = this->calcLength();
}

Edge::Edge(Node* n0, Node* n1, size_t value)
	: Element(MshElemType::LINE, value)
{
	_nodes = new Node*[2];
	_nodes[0] = n0;
	_nodes[1] = n1;

	this->_length = this->calcLength();
}

Edge::Edge(const Edge &edge)
	: Element(MshElemType::LINE, edge.getValue())
{
	_nodes = new Node*[2];
	_nodes[0] = edge._nodes[0];
	_nodes[1] = edge._nodes[1];
	_length = edge.getLength();
}

Edge::~Edge()
{
}

double Edge::calcLength()
{
	return sqrt(MathLib::sqrDist(_nodes[0]->getData(), _nodes[1]->getData()));
}

}

