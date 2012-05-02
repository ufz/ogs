/**
 * Quad.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Quad.h"
#include "Node.h"

#include "MathTools.h"
	
namespace MeshLib {

Quad::Quad(Node* nodes[4], size_t value)
	: Face(nodes, MshElemType::TRIANGLE, value)
{
	this->_area = this->calcArea();
}

Quad::Quad(Node* n0, Node* n1, Node* n2, Node* n3, size_t value)
	: Face(MshElemType::TRIANGLE, value)
{
	Node* nodes[4] = { n0, n1, n2, n3 };
	_nodes = nodes;

	this->_area = this->calcArea();
}

Quad::Quad(const Quad &quad)
	: Face(MshElemType::QUAD, quad.getValue())
{
	Node* nodes[4] = { new Node(*quad.getNode(0)), new Node(*quad.getNode(1)), new Node(*quad.getNode(2)), new Node(*quad.getNode(3)) };
	_area = quad.getArea();
}

Quad::~Quad()
{
}

double Quad::calcArea()
{
	return MathLib::calcDetTriangle(_nodes[0]->getData(), _nodes[1]->getData(), _nodes[2]->getData())
         + MathLib::calcDetTriangle(_nodes[2]->getData(), _nodes[3]->getData(), _nodes[0]->getData());
}

}