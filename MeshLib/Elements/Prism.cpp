/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file Prism.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Prism.h"
#include "Node.h"
#include "Tri.h"
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
	: Cell(MshElemType::PRISM, value)
{
	_nodes = nodes;
	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Prism::Prism(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, unsigned value)
	: Cell(MshElemType::PRISM, value)
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
	: Cell(MshElemType::PRISM, prism.getValue())
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

}

