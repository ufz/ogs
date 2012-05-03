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
	: Cell(nodes, MshElemType::PYRAMID, value)
{
	this->_volume = this->calcVolume();
}

Pyramid::Pyramid(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, size_t value)
	: Cell(MshElemType::PYRAMID, value)
{
	Node* nodes[5] = { n0, n1, n2, n3, n4 };
	_nodes = nodes;

	this->_volume = this->calcVolume();
}

Pyramid::Pyramid(const Pyramid &prism)
	: Cell(MshElemType::PYRAMID, prism.getValue())
{
	Node* nodes[5] = { new Node(*prism.getNode(0)), new Node(*prism.getNode(1)), new Node(*prism.getNode(2)), 
		               new Node(*prism.getNode(3)), new Node(*prism.getNode(4)) };
	_nodes = nodes;
	_volume = prism.getVolume();
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

