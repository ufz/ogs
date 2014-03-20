/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Element class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include "Element.h"
#include "Node.h"
#include "Line.h"

#include "MathTools.h"

namespace MeshLib {

Element::Element(unsigned value, std::size_t id)
	: _nodes(nullptr), _id(id), _value(value), _neighbors(nullptr)
{
}

Element::~Element()
{
	delete [] this->_nodes;
	delete [] this->_neighbors;
}

bool Element::addNeighbor(Element* e)
{
	if (e == this ||
		e == nullptr ||
		e->getDimension() != this->getDimension())
		return false;

	unsigned nNeighbors (this->getNNeighbors());
	for (unsigned n=0; n<nNeighbors; n++)
	{
		if (this->_neighbors[n] == e)
			return false;
		if (this->_neighbors[n] == nullptr)
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
		return new Line(nodes);
	}
	ERR("Error in MeshLib::Element::getEdge() - Index does not exist.");
	return nullptr;
}

void Element::computeSqrEdgeLengthRange(double &min, double &max) const
{
	min = std::numeric_limits<double>::max();
	max = 0;
	const unsigned nEdges (this->getNEdges());
	for (unsigned i=0; i<nEdges; i++)
	{
		const double dist (MathLib::sqrDist(*getEdgeNode(i,0), *getEdgeNode(i,1)));
		min = (dist<min) ? dist : min;
		max = (dist>max) ? dist : max;
	}
}

const Element* Element::getNeighbor(unsigned i) const
{
	if (i < getNNeighbors())
		return _neighbors[i];
	ERR("Error in MeshLib::Element::getNeighbor() - Index does not exist.");
	return nullptr;
}

unsigned Element::getNodeIDinElement(const MeshLib::Node* node) const
{
	const unsigned nNodes (this->getNNodes());
	for (unsigned i(0); i<nNodes; i++)
		if (node == _nodes[i])
			return i;
	return std::numeric_limits<unsigned>::max();
}

const Node* Element::getNode(unsigned i) const
{
	if (i < getNNodes())
		return _nodes[i];
	ERR("Error in MeshLib::Element::getNode() - Index %d in %s", i, MeshElemType2String(getGeomType()).c_str());
	return nullptr;
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
	ERR("Error in MeshLib::Element::getNodeIndex() - Index does not exist.");
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

