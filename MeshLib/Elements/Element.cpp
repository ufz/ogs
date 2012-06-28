/**
 * \file Element.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Element.h"
#include "Node.h"
#include "Edge.h"

#include "MathTools.h"

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

void Element::computeSqrEdgeLengthRange(double &min, double &max) const
{
	min = std::numeric_limits<double>::max();
	max = std::numeric_limits<double>::min();
	unsigned nEdges (this->getNEdges());
	for (unsigned i=0; i<nEdges; i++)
	{
		double dist (MathLib::sqrDist(getEdgeNode(i,0)->getCoords(), getEdgeNode(i,1)->getCoords()));
		min = (dist<min) ? dist : min;
		max = (dist>max) ? dist : max;
	}
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
	return std::numeric_limits<unsigned>::max();
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

