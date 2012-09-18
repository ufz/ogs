/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Prism.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Prism.h"
#include "Node.h"
#include "Tri.h"
#include "Pyramid.h"
#include "Quad.h"

#include "MathTools.h"

namespace MeshLib {

const unsigned Prism::_face_nodes[5][4] =
{
	{0, 2, 1, 99}, // Face 0
	{0, 1, 4,  3}, // Face 1
	{1, 2, 5,  4}, // Face 2
	{2, 0, 3,  5}, // Face 3
	{3, 4, 5, 99}  // Face 4
};

const unsigned Prism::_edge_nodes[9][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{0, 2}, // Edge 2
	{0, 3}, // Edge 3
	{1, 4}, // Edge 4
	{2, 5}, // Edge 5
	{3, 4}, // Edge 6
	{4, 5}, // Edge 7
	{3, 5}  // Edge 8
};

const unsigned Prism::_n_face_nodes[5] = { 3, 4, 4, 4, 3 };


Prism::Prism(Node* nodes[6], unsigned value)
	: Cell(value)
{
	_nodes = nodes;
	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Prism::Prism(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, unsigned value)
	: Cell(value)
{
	_nodes = new Node*[6];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_nodes[4] = n4;
	_nodes[5] = n5;
	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Prism::Prism(const Prism &prism)
	: Cell(prism.getValue())
{
	_nodes = new Node*[6];
	for (unsigned i=0; i<6; i++)
		_nodes[i] = prism._nodes[i];
	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++)
		_neighbors[i] = prism._neighbors[i];
	_volume = prism.getVolume();
}

Prism::~Prism()
{
}

double Prism::computeVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[3]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[1]->getCoords(), _nodes[4]->getCoords(), _nodes[2]->getCoords(), _nodes[3]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[2]->getCoords(), _nodes[4]->getCoords(), _nodes[5]->getCoords(), _nodes[3]->getCoords());
}

const Element* Prism::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];

		if (i==0 || i==4)
			return new Tri(nodes);
		else
			return new Quad(nodes);
	}
	std::cerr << "Error in MeshLib::Element::getFace() - Index does not exist." << std::endl;
	return NULL;
}

unsigned Prism::getNFaceNodes(unsigned i) const
{
	if (i<5)
		return _n_face_nodes[i];
	std::cerr << "Error in MeshLib::Element::getNFaceNodes() - Index does not exist." << std::endl;
	return 0;
}

bool Prism::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<9; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

Element* Prism::clone() const
{
	return new Prism(*this);
}

unsigned Prism::identifyFace(Node* nodes[3]) const
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

Element* Prism::reviseElement() const
{
	// try to create Pyramid
	if (_nodes[_edge_nodes[3][0]] == _nodes[_edge_nodes[3][1]]) {
		return new Pyramid(_nodes[1], _nodes[4], _nodes[5], _nodes[2], _nodes[0], _value);
	}

	if (_nodes[_edge_nodes[4][0]] == _nodes[_edge_nodes[4][1]]) {
		return new Pyramid(_nodes[0], _nodes[2], _nodes[5], _nodes[3], _nodes[1], _value);
	}

	if (_nodes[_edge_nodes[5][0]] == _nodes[_edge_nodes[5][1]]) {
		return new Pyramid(_nodes[0], _nodes[1], _nodes[4], _nodes[3], _nodes[2], _value);
	}

	return NULL;
}

} // end namespace MeshLib

