/**
 * Pyramid.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Pyramid.h"
#include "Node.h"

#include "MathTools.h"

namespace MeshLib {

Pyramid::Pyramid(Node* nodes[5], size_t value)
	: Cell(MshElemType::PYRAMID, value)
{
	_nodes = nodes;
	this->_volume = this->calcVolume();
}

Pyramid::Pyramid(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, size_t value)
	: Cell(MshElemType::PYRAMID, value)
{
	_nodes = new Node*[5];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_nodes[4] = n4;

	this->_volume = this->calcVolume();
}

Pyramid::Pyramid(const Pyramid &pyramid)
	: Cell(MshElemType::PYRAMID, pyramid.getValue())
{
	_nodes = new Node*[5];
	for (size_t i=0; i<5; i++)
		_nodes[i] = pyramid._nodes[i];
	_volume = pyramid.getVolume();
}

Pyramid::~Pyramid()
{
}

double Pyramid::calcVolume()
{
	return MathLib::calcDetTetrahedron(_nodes[0]->getData(), _nodes[1]->getData(), _nodes[2]->getData(), _nodes[4]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[2]->getData(), _nodes[3]->getData(), _nodes[0]->getData(), _nodes[4]->getData());
}

}

