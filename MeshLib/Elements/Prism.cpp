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
	: Cell(nodes, MshElemType::PRISM, value)
{
	this->_volume = this->calcVolume();
}

Prism::Prism(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, size_t value)
	: Cell(MshElemType::PRISM, value)
{
	Node* nodes[6] = { n0, n1, n2, n3, n4, n5 };
	_nodes = nodes;

	this->_volume = this->calcVolume();
}

Prism::Prism(const Prism &prism)
	: Cell(MshElemType::PRISM, prism.getValue())
{
	Node* nodes[6] = { new Node(*prism.getNode(0)), new Node(*prism.getNode(1)), new Node(*prism.getNode(2)), 
		               new Node(*prism.getNode(3)), new Node(*prism.getNode(4)), new Node(*prism.getNode(5)) };
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