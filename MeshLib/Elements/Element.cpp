/**
 * Element.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Element.h"
#include "Node.h"
#include "Edge.h"

#include <cassert>

namespace MeshLib {

/*
Element::Element(Node** nodes, MshElemType::type type, unsigned value)
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
}

const Element* Element::getEdge(unsigned i) const
{
	if (i < getNEdges())
	{
		Node** nodes = new Node*[2];
		nodes[0] = getEdgeNode(i,0);
		nodes[1] = getEdgeNode(i,1);
		return new Edge(nodes);
	}
	std::cerr << "Error in MeshLib::Element::getEdge() - Index does not exist." << std::endl;
	return NULL;
}

const Element* Element::getNeighbor(unsigned i) const
{
	if (i < getNNeighbors())
		return _neighbors[i];
	std::cerr << "Error in MeshLib::Element::getNeighbor() - Index does not exist." << std::endl;
	return NULL;
}

const Node* Element::getNode(unsigned i) const
{
	if (i < getNNodes())
		return _nodes[i];
	std::cerr << "Error in MeshLib::Element::getNode() - Index does not exist." << std::endl;
	return NULL;
}

unsigned Element::getNodeIndex(unsigned i) const 
{
	if (i<getNNodes())
		return _nodes[i]->getID();
	std::cerr << "Error in MeshLib::Element::getNodeIndex() - Index does not exist." << std::endl;
	return NULL;
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

