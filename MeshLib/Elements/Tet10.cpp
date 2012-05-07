/**
 * Tet10.cpp
 *
 *      Date: 2012/05/03
 *      Author: KR
 */

#include "Tet10.h"
#include "Node.h"

namespace MeshLib {

Tet10::Tet10(Node* nodes[10], size_t value)
	: Tet(value), FemElem()
{
	_nodes = nodes;
	this->_volume = this->calcVolume();
	this->calcCoG();
}

Tet10::Tet10(const Tet &tet)
	: Tet(tet.getValue()), FemElem()
{
	_nodes = new Node*[10];
	size_t nNodes (tet.getNNodes());
	for (size_t i=0; i<nNodes; i++)
		_nodes[i] = const_cast<Node*>(tet.getNode(i));

	if (nNodes < this->getNNodes())
	{
		//TODO: Calculate additional nodes!
	}

	this->_volume = this->calcVolume();
	this->calcCoG();
}

Tet10::Tet10(const Tet10 &tet)
	: Tet(tet.getValue()), FemElem()
{
	_nodes = new Node*[10];
	for (size_t i=0; i<10; i++)
		_nodes[i] = tet._nodes[i];
	_centre_of_gravity = tet.getCentreOfGravity();
}

Tet10::~Tet10()
{
}

void Tet10::calcCoG()
{
	//TODO calculation!
	_centre_of_gravity = 0;
}

}

