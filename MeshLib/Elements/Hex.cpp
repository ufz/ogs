/**
 * Hex.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Hex.h"
#include "Node.h"

#include "MathTools.h"


namespace MeshLib {

Hex::Hex(Node* nodes[8], size_t value)
	: Cell(nodes, MshElemType::HEXAHEDRON, value)
{
	this->_volume = this->calcVolume();
}

Hex::Hex(const Hex &hex)
	: Cell(MshElemType::HEXAHEDRON, hex.getValue())
{
	Node* nodes[8] = { new Node(*hex.getNode(0)), new Node(*hex.getNode(1)), new Node(*hex.getNode(2)), new Node(*hex.getNode(3)),
	                   new Node(*hex.getNode(4)), new Node(*hex.getNode(5)), new Node(*hex.getNode(6)), new Node(*hex.getNode(7)) };
	_nodes = nodes;
	_volume = hex.getVolume();
}

Hex::~Hex()
{
}

double Hex::calcVolume()
{
	return MathLib::calcDetTetrahedron(_nodes[4]->getData(), _nodes[7]->getData(), _nodes[5]->getData(), _nodes[0]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[5]->getData(), _nodes[3]->getData(), _nodes[1]->getData(), _nodes[0]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[5]->getData(), _nodes[7]->getData(), _nodes[3]->getData(), _nodes[0]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[5]->getData(), _nodes[7]->getData(), _nodes[6]->getData(), _nodes[2]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[1]->getData(), _nodes[3]->getData(), _nodes[5]->getData(), _nodes[2]->getData())
		 + MathLib::calcDetTetrahedron(_nodes[3]->getData(), _nodes[7]->getData(), _nodes[5]->getData(), _nodes[2]->getData());
}

}

