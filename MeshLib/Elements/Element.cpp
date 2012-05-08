/**
 * Element.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Element.h"
#include "Node.h"

#include <cassert>

namespace MeshLib {

/*
Element::Element(Node** nodes, MshElemType::type type, unsigned value)
	: _nodes(nodes), _type(type), _value(value)
{
}
*/
Element::Element(MshElemType::type type, unsigned value)
	: _nodes(NULL), _type(type), _value(value)
{
}

Element::~Element()
{
	delete[] this->_nodes;
	delete[] this->_neighbors;
}

const Element* Element::getNeighbor(unsigned i) const
{
	assert(i < getNNeighbors() && "Error in MeshLib::Element::getNeighbor() - Index does not exist.");
	return _neighbors[i];
}

const Node* Element::getNode(unsigned i) const
{
	assert(i < getNNodes() && "Error in MeshLib::Element::getNode() - Index does not exist.");
	assert(_nodes[i] != NULL && "Error in MeshLib::Element::getNode() - Node is NULL.");
	return _nodes[i];
}

unsigned Element::getNodeIndex(unsigned i) const 
{
	assert(i<getNNodes() && "Error in MeshLib::Element::getNodeIndex() - Index does not exist.");
	return _nodes[i]->getID();
}

bool Element::hasNeighbor(Element* elem) const
{
	unsigned nNeighbors (this->getNNeighbors());
	for (unsigned i=0; i<nNeighbors; i++)
		if (this->_neighbors[i]==elem)
			return true;
	return false;
}

}

