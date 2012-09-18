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
#include "Prism.h"

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
	: Cell(value)
{
	_nodes = nodes;
	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = NULL;
	this->_volume = this->computeVolume();
}

Hex::Hex(Node* n0, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, Node* n6, Node* n7, unsigned value)
	: Cell(value)
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
	: Cell(hex.getValue())
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

bool Hex::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<12; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

Element* Hex::clone() const
{
	return new Hex(*this);
}

unsigned Hex::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<6; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<4; j++)
			for (unsigned k=0; k<3; k++)
				if (_nodes[_face_nodes[i][j]] == nodes[k]) 
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

Element* Hex::reviseElement() const
{
	std::vector<size_t> collapsed_edges;
	for (size_t edge(0); edge<getNEdges(); edge++) {
		if (_nodes[_edge_nodes[edge][0]] == _nodes[_edge_nodes[edge][1]]) {
			collapsed_edges.push_back(edge);
		}
	}

	if (collapsed_edges.size() == 1) {
		std::cerr << "[Hex::reviseElement()] collapsing of one edge in hexahedron not handled" << std::endl;
		return NULL;
	}

	if (collapsed_edges.size() == 2) {
		// try to create a prism out of the hex
		if (collapsed_edges[0] == 0 && collapsed_edges[1] == 2) {
			return new Prism(_nodes[0],_nodes[4], _nodes[5], _nodes[3], _nodes[7], _nodes[6], _value);
		}
		if (collapsed_edges[0] == 1 && collapsed_edges[1] == 3) {
			return new Prism(_nodes[0],_nodes[4], _nodes[7], _nodes[1], _nodes[5], _nodes[6], _value);
		}
		if (collapsed_edges[0] == 4 && collapsed_edges[1] == 5) {
			return new Prism(_nodes[0],_nodes[7], _nodes[3], _nodes[1], _nodes[6], _nodes[2], _value);
		}
		if (collapsed_edges[0] == 5 && collapsed_edges[1] == 6) {
			return new Prism(_nodes[0],_nodes[1], _nodes[4], _nodes[3], _nodes[2], _nodes[7], _value);
		}
		if (collapsed_edges[0] == 6 && collapsed_edges[1] == 7) {
			return new Prism(_nodes[0],_nodes[3], _nodes[4], _nodes[1], _nodes[2], _nodes[5], _value);
		}
		if (collapsed_edges[0] == 7 && collapsed_edges[1] == 4) {
			return new Prism(_nodes[0],_nodes[1], _nodes[5], _nodes[3], _nodes[2], _nodes[6], _value);
		}
		if (collapsed_edges[0] == 8 && collapsed_edges[1] == 10) {
			return new Prism(_nodes[0],_nodes[1], _nodes[4], _nodes[3], _nodes[2], _nodes[7], _value);
		}
		if (collapsed_edges[0] == 9 && collapsed_edges[1] == 11) {
			return new Prism(_nodes[0],_nodes[3], _nodes[4], _nodes[1], _nodes[2], _nodes[5], _value);
		}
		return NULL;
	}

	return NULL;
}

} // end namespace MeshLib
