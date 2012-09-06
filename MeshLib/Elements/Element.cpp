/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
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
Element::Element(unsigned value)
	: _nodes(NULL), _value(value), _neighbors(NULL)
{
}

Element::~Element()
{
	delete [] this->_nodes;
	delete [] this->_neighbors;
}

bool Element::addNeighbor(Element* e)
{
	if (e == this)
		return false;

	unsigned nNeighbors (this->getNNeighbors());
	for (unsigned n=0; n<nNeighbors; n++)
	{
		if (this->_neighbors[n] == e)
			return false;
		if (this->_neighbors[n] == NULL)
			break;
	}

	Node* face_nodes[3];
	const unsigned nNodes (this->getNNodes());
	const unsigned eNodes (e->getNNodes());
	const Node* const* e_nodes = e->getNodes();
	unsigned count(0);
	const unsigned dim (this->getDimension());
	for (unsigned i(0); i<nNodes; i++)
		for (unsigned j(0); j<eNodes; j++)
			if (_nodes[i] == e_nodes[j])
			{
				face_nodes[count] = _nodes[i];
				//std::cout << _nodes[i]->getID() << " == " << e_nodes[j]->getID() << std::endl;
				// increment shared nodes counter and check if enough nodes are similar to be sure e is a neighbour of this
				if ((++count)>=dim)
				{
					_neighbors[ this->identifyFace(face_nodes) ] = e;
					return true;
				}
			}

	return false;
}

const Element* Element::getEdge(unsigned i) const
{
	if (i < getNEdges())
	{
		Node** nodes = new Node*[2];
		nodes[0] = const_cast<Node*>(getEdgeNode(i,0));
		nodes[1] = const_cast<Node*>(getEdgeNode(i,1));
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
		double dist (MathLib::sqrDist(getEdgeNode(i,0), getEdgeNode(i,1)));
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
	std::cerr << "Error in MeshLib::Element::getNode() - Index " << i << " in " << MshElemType2String(getType()) <<  " does not exist." << std::endl;
	return NULL;
}

void Element::setNode(unsigned idx, Node* node)
{
	if (idx < getNNodes())
		_nodes[idx] = node;
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

