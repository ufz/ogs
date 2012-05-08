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
	: Face(MshElemType::TRIANGLE, value)
{
	_nodes = nodes;
	_neighbors = new Element*[4];
	for (size_t i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_area = this->calcArea();
}

Quad::Quad(Node* n0, Node* n1, Node* n2, Node* n3, size_t value)
	: Face(MshElemType::TRIANGLE, value)
{
	_nodes = new Node*[4];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_neighbors = new Element*[4];
	for (size_t i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_area = this->calcArea();
}

Quad::Quad(const Quad &quad)
	: Face(MshElemType::QUAD, quad.getValue())
{
	_nodes = new Node*[4];
	_neighbors = new Element*[4];
	for (size_t i=0; i<4; i++)
	{
		_nodes[i] = quad._nodes[i];
		_neighbors[i] = quad._neighbors[i];
	}
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

