/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tet.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Tet.h"
#include "Node.h"
#include "Tri.h"

#include "MathTools.h"

namespace MeshLib {


const unsigned Tet::_face_nodes[4][3] =
{
	{0, 2, 1}, // Face 0
	{0, 1, 3}, // Face 1
	{1, 2, 3}, // Face 2
	{2, 0, 3}  // Face 3
};

const unsigned Tet::_edge_nodes[6][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{0, 2}, // Edge 2
	{0, 3}, // Edge 3
	{1, 3}, // Edge 4
	{2, 3}  // Edge 5
};


Tet::Tet(Node* nodes[4], unsigned value)
	: Cell(value)
{
	_nodes = nodes;
	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Tet::Tet(Node* n0, Node* n1, Node* n2, Node* n3, unsigned value)
	: Cell(value)
{
	_nodes = new Node*[4];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Tet::Tet(unsigned value)
	: Cell(value)
{
	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
		_neighbors[i] = NULL;
}

Tet::Tet(const Tet &tet)
	: Cell(tet.getValue())
{
	_nodes = new Node*[4];
	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
	{
		_nodes[i] = tet._nodes[i];
		_neighbors[i] = tet._neighbors[i];
	}
	_volume = tet.getVolume();
}

Tet::~Tet()
{
}

double Tet::computeVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[3]->getCoords());
}

const Element* Tet::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];
		return new Tri(nodes);
	}
	std::cerr << "Error in MeshLib::Element::getFace() - Index does not exist." << std::endl;
	return NULL;
}

Element* Tet::clone() const
{
	return new Tet(*this);
}


}

