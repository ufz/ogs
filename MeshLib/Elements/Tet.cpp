/**
 * Tet.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Tet.h"
#include "Node.h"

#include "MathTools.h"

namespace MeshLib {

Tet::Tet(Node* nodes[4], size_t value)
	: Cell(nodes, MshElemType::TETRAHEDRON, value)
{
	this->_volume = this->calcVolume();
}

Tet::Tet(Node* n0, Node* n1, Node* n2, Node* n3, size_t value)
	: Cell(MshElemType::TETRAHEDRON, value)
{
	Node* nodes[4] = { n0, n1, n2, n3 };
	_nodes = nodes;

	this->_volume = this->calcVolume();
}

Tet::Tet(size_t value)
	: Cell(MshElemType::TETRAHEDRON, value)
{
}

Tet::Tet(const Tet &tet)
	: Cell(MshElemType::TETRAHEDRON, tet.getValue())
{
	Node* nodes[4] = { new Node(*tet.getNode(0)), new Node(*tet.getNode(1)), new Node(*tet.getNode(2)), new Node(*tet.getNode(3)) };
	_nodes = nodes;
	_volume = tet.getVolume();
}

Tet::~Tet()
{
}

double Tet::calcVolume()
{
	return MathLib::calcDetTetrahedron(_nodes[0]->getData(), _nodes[1]->getData(), _nodes[2]->getData(), _nodes[3]->getData());
}

}

