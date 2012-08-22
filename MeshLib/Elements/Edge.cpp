/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Edge.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Edge.h"
#include "Node.h"

#include "MathTools.h"

namespace MeshLib {

Edge::Edge(Node* nodes[2], unsigned value)
	: Element(value)
{
	_nodes = nodes;
	this->_length = this->computeLength();
}

Edge::Edge(Node* n0, Node* n1, unsigned value)
	: Element(value)
{
	_nodes = new Node*[2];
	_nodes[0] = n0;
	_nodes[1] = n1;

	this->_length = this->computeLength();
}

Edge::Edge(const Edge &edge)
	: Element(edge.getValue())
{
	_nodes = new Node*[2];
	_nodes[0] = edge._nodes[0];
	_nodes[1] = edge._nodes[1];
	_length = edge.getLength();
}

Edge::~Edge()
{
}

double Edge::computeLength()
{
	return sqrt(MathLib::sqrDist(_nodes[0]->getCoords(), _nodes[1]->getCoords()));
}

Element* Edge::clone() const
{
	return new Edge(*this);
}

Element* Edge::reviseElement() const
{
	if (_nodes[0] == _nodes[1]) {
		return NULL;
	}

	return NULL;
}

}

