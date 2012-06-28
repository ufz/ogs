/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Hex.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Hex.h"
#include "Node.h"
#include "Quad.h"

#include "MathTools.h"


namespace MeshLib {

const unsigned Hex::_face_nodes[6][4] =
{
	{0, 3, 2, 1}, // Face 0
	{0, 1, 5, 4}, // Face 1
	{1, 2, 6, 5}, // Face 2
	{2, 3, 7, 6}, // Face 3
	{3, 0, 4, 7}, // Face 4
	{4, 5, 6, 7}  // Face 5
};

const unsigned Hex::_edge_nodes[12][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}, // Edge 3
	{0, 4}, // Edge 4
	{1, 5}, // Edge 5
	{2, 6}, // Edge 6
	{3, 7}, // Edge 7
	{4, 5}, // Edge 8
	{5, 6}, // Edge 9
	{6, 7}, // Edge 10
	{4, 7}  // Edge 11
};


Hex::Hex(Node* nodes[8], unsigned value)
	: Cell(MshElemType::HEXAHEDRON, value)
{
	_nodes = nodes;
	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Hex::Hex(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, Node* n6, Node* n7, unsigned value)
	: Cell(MshElemType::HEXAHEDRON, value)
{
	_nodes = new Node*[8];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_nodes[3] = n3;
	_nodes[4] = n4;
	_nodes[5] = n5;
	_nodes[6] = n6;
	_nodes[7] = n7;
	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Hex::Hex(const Hex &hex)
	: Cell(MshElemType::HEXAHEDRON, hex.getValue())
{
	_nodes = new Node*[8];
	for (unsigned i=0; i<8; i++)
		_nodes[i] = hex._nodes[i];
	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = hex._neighbors[i];
	_volume = hex.getVolume();
}

Hex::~Hex()
{
}

double Hex::computeVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[4]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[0]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[3]->getCoords(), _nodes[1]->getCoords(), _nodes[0]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[6]->getCoords(), _nodes[2]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[1]->getCoords(), _nodes[3]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[3]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords());
}

const Element* Hex::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];
		return new Quad(nodes);
	}
	std::cerr << "Error in MeshLib::Element::getFace() - Index does not exist." << std::endl;
	return NULL;
}

}

