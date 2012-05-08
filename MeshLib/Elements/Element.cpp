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
Element::Element(Node** nodes, MshElemType::type type, size_t value)
	: _nodes(nodes), _type(type), _value(value)
{
}
*/
Element::Element(MshElemType::type type, size_t value)
	: _nodes(NULL), _type(type), _value(value)
{
}

Element::~Element()
{
	delete[] this->_nodes;
	delete[] this->_neighbors;
}

const Element* Element::getNeighbor(size_t i) const
{
	assert(i < getNNeighbors() && "Error in MeshLib::Element::getNeighbor() - Index does not exist.");
	return _neighbors[i];
}

const Node* Element::getNode(size_t i) const
{
	assert(i < getNNodes() && "Error in MeshLib::Element::getNode() - Index does not exist.");
	assert(_nodes[i] != NULL && "Error in MeshLib::Element::getNode() - Node is NULL.");
	return _nodes[i];
}

size_t Element::getNodeIndex(size_t i) const 
{
	assert(i<getNNodes() && "Error in MeshLib::Element::getNodeIndex() - Index does not exist.");
	return _nodes[i]->getID();
}

bool Element::hasNeighbor(Element* elem) const
{
	size_t nNeighbors (this->getNNeighbors());
	for (size_t i=0; i<nNeighbors; i++)
		if (this->_neighbors[i]==elem)
			return true;
	return false;
}

}

