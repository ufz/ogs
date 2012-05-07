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
	: Cell(MshElemType::HEXAHEDRON, value)
{
	_nodes = nodes;
	this->_volume = this->calcVolume();
}

Hex::Hex(const Hex &hex)
	: Cell(MshElemType::HEXAHEDRON, hex.getValue())
{
	_nodes = new Node*[8];
	for (size_t i=0; i<8; i++)
		_nodes[i] = hex._nodes[i];
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

