/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Pyramid.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Pyramid.h"
#include "Node.h"
#include "Tri.h"
#include "Tet.h"
#include "Quad.h"

#include "MathTools.h"

namespace MeshLib {

const unsigned Pyramid::_face_nodes[5][4] =
{
	{0, 1, 4, 99}, // Face 0
	{1, 2, 4, 99}, // Face 1
	{2, 3, 4, 99}, // Face 2
	{3, 0, 4, 99}, // Face 3
	{0, 3, 2,  1}  // Face 4
};

const unsigned Pyramid::_edge_nodes[8][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}, // Edge 3
	{0, 4}, // Edge 4
	{1, 4}, // Edge 5
	{2, 4}, // Edge 6
	{3, 4}  // Edge 7
};

const unsigned Pyramid::_n_face_nodes[5] = { 3, 3, 3, 3, 4 };


Pyramid::Pyramid(Node* nodes[5], unsigned value)
	: Cell(value)
{
	_nodes = nodes;
	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Pyramid::Pyramid(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, unsigned value)
	: Cell(value)
{
	_nodes = new Node*[5];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_nodes[4] = n4;
	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++)
		_neighbors[i] = NULL;

	this->_volume = this->computeVolume();
}

Pyramid::Pyramid(const Pyramid &pyramid)
	: Cell(pyramid.getValue())
{
	_nodes = new Node*[5];
	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++)
	{
		_nodes[i] = pyramid._nodes[i];
		_neighbors[i] = pyramid._neighbors[i];
	}
	_volume = pyramid.getVolume();
}

Pyramid::~Pyramid()
{
}

double Pyramid::computeVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[4]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[2]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords(), _nodes[4]->getCoords());
}

const Element* Pyramid::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];

		if (i<4)
			return new Tri(nodes);
		else
			return new Quad(nodes);
	}
	std::cerr << "Error in MeshLib::Element::getFace() - Index does not exist." << std::endl;
	return NULL;
}

unsigned Pyramid::getNFaceNodes(unsigned i) const
{
	if (i<5)
		return _n_face_nodes[i];
	std::cerr << "Error in MeshLib::Element::getNFaceNodes() - Index does not exist." << std::endl;
	return 0;
}

bool Pyramid::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<8; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

Element* Pyramid::clone() const
{
	return new Pyramid(*this);
}

unsigned Pyramid::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<5; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<4; j++)
			for (unsigned k=0; k<3; k++)
				if (_face_nodes[i][j] != 99 && _nodes[_face_nodes[i][j]] == nodes[k]) 
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

Element* Pyramid::reviseElement() const
{
	// try to create tetrahedron
	if (_nodes[_edge_nodes[0][0]] == _nodes[_edge_nodes[0][1]]
		|| _nodes[_edge_nodes[1][0]] == _nodes[_edge_nodes[1][1]]) {
		return new Tet(_nodes[0], _nodes[2], _nodes[3], _nodes[4], _value);
	}

	if (_nodes[_edge_nodes[2][0]] == _nodes[_edge_nodes[2][1]]
	|| _nodes[_edge_nodes[3][0]] == _nodes[_edge_nodes[3][1]]) {
		return new Tet(_nodes[0], _nodes[1], _nodes[2], _nodes[4], _value);
	}

	return NULL;
}

} // end namespace MeshLib
