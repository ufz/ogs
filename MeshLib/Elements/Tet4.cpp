/**
 * Tet4.cpp
 *
 *      Date: 2012/05/03
 *      Author: KR
 */

#include "Tet4.h"
#include "Node.h"

namespace MeshLib {

Tet4::Tet4(Node* nodes[4], size_t value)
	: Tet(nodes, value), FemElem()
{
	this->calcCoG();
}

Tet4::Tet4(const Tet &tet)
	: Tet(tet), FemElem()
{
	this->calcCoG();
}

Tet4::Tet4(const Tet4 &tet)
	: Tet(tet.getValue()), FemElem()
{
	for (size_t i=0; i<4; i++)
		_nodes[i] = tet._nodes[i];
	_centre_of_gravity = tet.getCentreOfGravity();
}

Tet4::~Tet4()
{
}

void Tet4::calcCoG()
{
	_centre_of_gravity = 0;
	//TODO calculation!
}

}

