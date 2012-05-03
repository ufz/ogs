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
	_centre_of_gravity = this->calcCoG();
}

Tet10::Tet10(const Tet &tet)
	: Tet(tet.getValue()), FemElem()
{
	Node* nodes[10]; //= { n0, n1, n2, n3 };
	//TODO: Calculate additional nodes!
	_nodes = nodes;
	this->_volume = this->calcVolume();
	_centre_of_gravity = this->calcCoG();
}

Tet10::Tet10(const Tet10 &tet)
	: Tet(tet.getValue()), FemElem()
{
	Node* nodes[10] = { new Node(*tet.getNode(0)), new Node(*tet.getNode(1)), new Node(*tet.getNode(2)), 
		                new Node(*tet.getNode(3)), new Node(*tet.getNode(4)), new Node(*tet.getNode(5)), 
						new Node(*tet.getNode(6)), new Node(*tet.getNode(7)), new Node(*tet.getNode(8)), 
						new Node(*tet.getNode(9)) };
	_nodes = nodes;
	_centre_of_gravity = tet.getCentreOfGravity();
}

Tet10::~Tet10()
{
}

GEOLIB::Point Tet10::calcCoG()
{
	GEOLIB::Point cog(0,0,0);
	//TODO calculation!
	return cog;
}

}

