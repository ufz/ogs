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
	_centre_of_gravity = this->calcCoG();
}

Tet4::Tet4(const Tet &tet)
	: Tet(tet.getValue()), FemElem()
{
	Node* nodes[4] = { new Node(*tet.getNode(0)), new Node(*tet.getNode(1)), new Node(*tet.getNode(2)), new Node(*tet.getNode(3)) };
	_nodes = nodes;
	_centre_of_gravity = this->calcCoG();
}

Tet4::Tet4(const Tet4 &tet)
	: Tet(tet.getValue()), FemElem()
{
	Node* nodes[4] = { new Node(*tet.getNode(0)), new Node(*tet.getNode(1)), new Node(*tet.getNode(2)), new Node(*tet.getNode(3)) };
	_nodes = nodes;
	_centre_of_gravity = tet.getCentreOfGravity();
}

Tet4::~Tet4()
{
}

GEOLIB::Point Tet4::calcCoG()
{
	GEOLIB::Point cog(0,0,0);
	//TODO calculation!
	return cog;
}

}

