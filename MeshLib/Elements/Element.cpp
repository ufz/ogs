/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Element class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

void Element::setNeighbor(Element* neighbor, unsigned const face_id)
{
	if (neighbor == this)
		return;

	this->_neighbors[face_id] = neighbor;
}

boost::optional<unsigned> Element::addNeighbor(Element* e)
{
	if (e == this ||
		e == nullptr ||
		e->getDimension() != this->getDimension())
		return boost::optional<unsigned>();

	if (this->hasNeighbor(e))
		return boost::optional<unsigned>();

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
					return boost::optional<unsigned>(e->identifyFace(face_nodes));
				}
			}

	return boost::optional<unsigned>();
}

MeshLib::Node Element::getCenterOfGravity() const
{
	const unsigned nNodes (this->getNNodes());
	MeshLib::Node center(0,0,0);
	for (unsigned i=0; i<nNodes; ++i)
	{
		center[0] += (*_nodes[i])[0];
		center[1] += (*_nodes[i])[1];
		center[2] += (*_nodes[i])[2];
	}
	center[0] /= nNodes;
	center[1] /= nNodes;
	center[2] /= nNodes;
	return center;
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

bool Element::isBoundaryElement() const
{
    return std::any_of(_neighbors, _neighbors + this->getNNeighbors(), 
        [](MeshLib::Element* e){ return e == nullptr; });
}

}

