/**
 * Tri.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Tri.h"
#include "Node.h"

#include "MathTools.h"

namespace MeshLib {

Tri::Tri(Node* nodes[3], size_t value)
	: Face(nodes, MshElemType::TRIANGLE, value)
{
	this->_area = this->calcArea();
}

Tri::Tri(Node* n0, Node* n1, Node* n2, size_t value)
	: Face(MshElemType::TRIANGLE, value)
{
	Node* nodes[3] = { n0, n1, n2 };
	_nodes = nodes;

	this->_area = this->calcArea();
}

Tri::Tri(const Tri &tri)
	: Face(MshElemType::TRIANGLE, tri.getValue())
{
	Node* nodes[3] = { new Node(*tri.getNode(0)), new Node(*tri.getNode(1)), new Node(*tri.getNode(2)) };
	_area = tri.getArea();
}

Tri::~Tri()
{
}

double Tri::calcArea()
{
	return MathLib::calcDetTriangle(_nodes[0]->getData(), _nodes[1]->getData(), _nodes[2]->getData());
}

}