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
	: Face(MshElemType::TRIANGLE, value)
{
	_nodes = nodes;
	_neighbors = new Element*[3];
	for (size_t i=0; i<3; i++)
		_neighbors[i] = NULL;
	this->_area = this->calcArea();
}

Tri::Tri(Node* n0, Node* n1, Node* n2, size_t value)
	: Face(MshElemType::TRIANGLE, value)
{
	_nodes = new Node*[3];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_neighbors = new Element*[3];
	for (size_t i=0; i<3; i++)
		_neighbors[i] = NULL;
	this->_area = this->calcArea();
}

Tri::Tri(const Tri &tri)
	: Face(MshElemType::TRIANGLE, tri.getValue())
{
	_nodes = new Node*[3];
	_neighbors = new Element*[3];
	for (size_t i=0; i<3; i++)
	{
		_nodes[i] = tri._nodes[i];
		_neighbors[i] = tri._neighbors[i];
	}
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

