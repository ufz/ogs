/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Quad.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Quad.h"
#include "Node.h"
#include "Tri.h"

#include "MathTools.h"

namespace MeshLib {


const unsigned Quad::_edge_nodes[4][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}  // Edge 3
};


Quad::Quad(Node* nodes[4], unsigned value)
	: Face(value)
{
	_nodes = nodes;
	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_area = this->computeVolume();
}

Quad::Quad(Node* n0, Node* n1, Node* n2, Node* n3, unsigned value)
	: Face(value)
{
	_nodes = new Node*[4];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_area = this->computeVolume();
}

Quad::Quad(const Quad &quad)
	: Face(quad.getValue())
{
	_nodes = new Node*[4];
	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
	{
		_nodes[i] = quad._nodes[i];
		_neighbors[i] = quad._neighbors[i];
	}
	_area = quad.getArea();
}

Quad::~Quad()
{
}

double Quad::computeVolume()
{
	return MathLib::calcTriangleArea(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords())
         + MathLib::calcTriangleArea(_nodes[2]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords());
}

Element* Quad::clone() const
{
	return new Quad(*this);
}

unsigned Quad::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<4; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<2; j++)
			for (unsigned k=0; k<2; k++)
				if (_nodes[_edge_nodes[i][j]] == nodes[k]) 
					flag++;
		if (flag==2)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}
Element* Quad::reviseElement() const
{
	if (_nodes[0] == _nodes[1] || _nodes[1] == _nodes[2]) {
		return new Tri(_nodes[0], _nodes[2], _nodes[3], _value);
	}

	if (_nodes[2] == _nodes[3] || _nodes[3] == _nodes[0]) {
		return new Tri(_nodes[0], _nodes[1], _nodes[2], _value);
	}

	// this should not happen
	return NULL;
}

}

