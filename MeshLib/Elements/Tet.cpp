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
	: Cell(MshElemType::TETRAHEDRON, value)
{
	_nodes = nodes;
	_neighbors = new Element*[4];
	for (size_t i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_volume = this->calcVolume();
}

Tet::Tet(Node* n0, Node* n1, Node* n2, Node* n3, size_t value)
	: Cell(MshElemType::TETRAHEDRON, value)
{
	_nodes = new Node*[4];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_neighbors = new Element*[4];
	for (size_t i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_volume = this->calcVolume();
}

Tet::Tet(size_t value)
	: Cell(MshElemType::TETRAHEDRON, value)
{
	_neighbors = new Element*[4];
	for (size_t i=0; i<4; i++)
		_neighbors[i] = NULL;
}

Tet::Tet(const Tet &tet)
	: Cell(MshElemType::TETRAHEDRON, tet.getValue())
{
	_nodes = new Node*[4];
	_neighbors = new Element*[4];
	for (size_t i=0; i<4; i++)
	{
		_nodes[i] = tet._nodes[i];
		_neighbors[i] = tet._neighbors[i];
	}
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

