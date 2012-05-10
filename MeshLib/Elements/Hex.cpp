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

Hex::Hex(Node* nodes[8], unsigned value)
	: Cell(MshElemType::HEXAHEDRON, value)
{
	_nodes = nodes;
	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = NULL;
	this->_volume = this->calcVolume();
}

Hex::Hex(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, Node* n6, Node* n7, unsigned value)
	: Cell(MshElemType::HEXAHEDRON, value)
{
	_nodes = new Node*[8];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_nodes[4] = n4;
	_nodes[5] = n5;
	_nodes[6] = n6;
	_nodes[7] = n7;
	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = NULL;
	this->_volume = this->calcVolume();
}

Hex::Hex(const Hex &hex)
	: Cell(MshElemType::HEXAHEDRON, hex.getValue())
{
	_nodes = new Node*[8];
	for (unsigned i=0; i<8; i++)
		_nodes[i] = hex._nodes[i];
	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = hex._neighbors[i];
	_volume = hex.getVolume();
}

Hex::~Hex()
{
}

double Hex::calcVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[4]->getData(), _nodes[7]->getData(), _nodes[5]->getData(), _nodes[0]->getData())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getData(), _nodes[3]->getData(), _nodes[1]->getData(), _nodes[0]->getData())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getData(), _nodes[7]->getData(), _nodes[3]->getData(), _nodes[0]->getData())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getData(), _nodes[7]->getData(), _nodes[6]->getData(), _nodes[2]->getData())
		 + MathLib::calcTetrahedronVolume(_nodes[1]->getData(), _nodes[3]->getData(), _nodes[5]->getData(), _nodes[2]->getData())
		 + MathLib::calcTetrahedronVolume(_nodes[3]->getData(), _nodes[7]->getData(), _nodes[5]->getData(), _nodes[2]->getData());
}

}

