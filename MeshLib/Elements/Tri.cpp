/**
 * \file Tri.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Tri.h"
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
	: Face(MshElemType::TRIANGLE, value)
{
	_nodes = nodes;
	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++)
		_neighbors[i] = NULL;
	this->_area = this->computeArea();
}

Tri::Tri(Node* n0, Node* n1, Node* n2, unsigned value)
	: Face(MshElemType::TRIANGLE, value)
{
	_nodes = new Node*[3];
	_nodes[0] = n0;
	_nodes[1] = n1;
	_nodes[2] = n2;
	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++)
		_neighbors[i] = NULL;
	this->_area = this->computeArea();
}

Tri::Tri(const Tri &tri)
	: Face(MshElemType::TRIANGLE, tri.getValue())
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

double Tri::computeArea()
{
	return MathLib::calcTriangleArea(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords());
}

}

