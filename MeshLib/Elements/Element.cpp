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

Element::Element(Node** nodes, MshElemType::type type, size_t value)
	: _nodes(nodes), _type(type), _value(value)
{
}

Element::Element(MshElemType::type type, size_t value)
	: _type(type), _value(value)
{
}

Element::~Element()
{
	delete[] this->_nodes;
}

const Node* Element::getNode(size_t i) const
{
	assert(i<getNNodes() && "Error in MeshLib::Element - Index does not exist.");
	return _nodes[i];
}

size_t Element::getNodeIndex(size_t i) const 
{
	assert(i<getNNodes() && "Error in MeshLib::Element - Index does not exist.");
	return _nodes[i]->getID();
}

}