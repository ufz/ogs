/**
 * Prism.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Prism.h"
#include "Node.h"

#include "MathTools.h"

namespace MeshLib {

Prism::Prism(Node* nodes[6], size_t value)
	: Cell(MshElemType::PRISM, value)
{
	_nodes = _nodes;
	_neighbors = new Element*[5];
	for (size_t i=0; i<5; i++)
		_neighbors[i] = NULL;
	this->_volume = this->calcVolume();
}

Prism::Prism(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, size_t value)
	: Cell(MshElemType::PRISM, value)
{
	_nodes = new Node*[6];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_nodes[4] = n4;
	_nodes[5] = n5;
	_neighbors = new Element*[5];
	for (size_t i=0; i<5; i++)
		_neighbors[i] = NULL;
	this->_volume = this->calcVolume();
}

Prism::Prism(const Prism &prism)
	: Cell(MshElemType::PRISM, prism.getValue())
{
	_nodes = new Node*[6];
	for (size_t i=0; i<6; i++)
		_nodes[i] = prism._nodes[i];
	_neighbors = new Element*[5];
	for (size_t i=0; i<5; i++)
		_neighbors[i] = prism._neighbors[i];
	_volume = prism.getVolume();
}

Prism::~Prism()
{
}

double Prism::calcVolume()
{
	return MathLib::calcDetTetrahedron(_nodes[0]->getData(), _nodes[1]->getData(), _nodes[2]->getData(), _nodes[3]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[1]->getData(), _nodes[4]->getData(), _nodes[2]->getData(), _nodes[3]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[2]->getData(), _nodes[4]->getData(), _nodes[5]->getData(), _nodes[3]->getData());
}

}

