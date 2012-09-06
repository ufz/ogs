/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tri.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Tri.h"
#include "Edge.h"
#include "Node.h"

#include "MathTools.h"

namespace MeshLib {


const unsigned Tri::_edge_nodes[3][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{0, 2}  // Edge 2
};


Tri::Tri(Node* nodes[3], unsigned value)
	: Face(value)
{
	_nodes = nodes;
	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++)
		_neighbors[i] = NULL;
	this->_area = this->computeVolume();
}

Tri::Tri(Node* n0, Node* n1, Node* n2, unsigned value)
	: Face(value)
{
	_nodes = new Node*[3];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++)
		_neighbors[i] = NULL;
	this->_area = this->computeVolume();
}

Tri::Tri(const Tri &tri)
	: Face(tri.getValue())
{
	_nodes = new Node*[3];
	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++)
	{
		_nodes[i] = tri._nodes[i];
		_neighbors[i] = tri._neighbors[i];
	}
	_area = tri.getArea();
}

Tri::~Tri()
{
}

Element* Tri::clone() const
{
	return new Tri(*this);
}

double Tri::computeVolume()
{
	return MathLib::calcTriangleArea(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords());
}

unsigned Tri::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<3; i++)
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

Element* Tri::reviseElement() const
{
	// try to create an edge
	if (_nodes[0] == _nodes[1] || _nodes[1] == _nodes[2]) {
		return new Edge(_nodes[0], _nodes[2], _value);
	}

	if (_nodes[0] == _nodes[2]) {
		return new Edge(_nodes[0], _nodes[1], _value);
	}

	return NULL;
}

}

