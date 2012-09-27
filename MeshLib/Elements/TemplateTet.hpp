/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file TemplateTet.hpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Node.h"
#include "Tri.h"

#include "MathTools.h"

namespace MeshLib {

template <unsigned ORDER, unsigned NNODES>
const unsigned TemplateTet<ORDER,NNODES>::_face_nodes[4][3] =
{
	{0, 2, 1}, // Face 0
	{0, 1, 3}, // Face 1
	{1, 2, 3}, // Face 2
	{2, 0, 3}  // Face 3
};

template <unsigned ORDER, unsigned NNODES>
const unsigned TemplateTet<ORDER, NNODES>::_edge_nodes[6][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{0, 2}, // Edge 2
	{0, 3}, // Edge 3
	{1, 3}, // Edge 4
	{2, 3}  // Edge 5
};

template <unsigned ORDER, unsigned NNODES>
TemplateTet<ORDER, NNODES>::TemplateTet(Node* nodes[NNODES], unsigned value)
	: Cell(value)
{
	_nodes = nodes;

	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
		_neighbors[i] = NULL;

	this->_volume = this->computeVolume();
}

template <unsigned ORDER, unsigned NNODES>
TemplateTet<ORDER, NNODES>::TemplateTet(const TemplateTet<ORDER, NNODES> &tet)
	: Cell(tet.getValue())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++) {
		_nodes[i] = tet._nodes[i];
	}

	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
	{
		_neighbors[i] = tet._neighbors[i];
	}

	_volume = tet.getVolume();
}

template <unsigned ORDER, unsigned NNODES>
TemplateTet<ORDER, NNODES>::~TemplateTet()
{
}

template <unsigned ORDER, unsigned NNODES>
double TemplateTet<ORDER, NNODES>::computeVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[3]->getCoords());
}

template <unsigned ORDER, unsigned NNODES>
const Element* TemplateTet<ORDER, NNODES>::getFace(unsigned i) const
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

template <unsigned ORDER, unsigned NNODES>
bool TemplateTet<ORDER, NNODES>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<6; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned ORDER, unsigned NNODES>
Element* TemplateTet<ORDER, NNODES>::clone() const
{
	return new TemplateTet<ORDER,NNODES>(*this);
}

template <unsigned ORDER, unsigned NNODES>
unsigned TemplateTet<ORDER, NNODES>::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<4; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<3; j++)
			for (unsigned k=0; k<3; k++)
				if (_nodes[_face_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

template <unsigned ORDER, unsigned NNODES>
Element* TemplateTet<ORDER, NNODES>::reviseElement() const
{
	if (_nodes[0] == _nodes[1] || _nodes[1] == _nodes[2]) {
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		tri_nodes[0] = _nodes[0];
		tri_nodes[1] = _nodes[2];
		tri_nodes[2] = _nodes[3];
		return new Tri(tri_nodes, _value);
	}

	if (_nodes[2] == _nodes[0]) {
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		tri_nodes[0] = _nodes[0];
		tri_nodes[1] = _nodes[1];
		tri_nodes[2] = _nodes[3];
		return new Tri(tri_nodes, _value);
	}

	if (_nodes[0] == _nodes[3] || _nodes[1] == _nodes[3] || _nodes[2] == _nodes[3]) {
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		tri_nodes[0] = _nodes[0];
		tri_nodes[1] = _nodes[1];
		tri_nodes[2] = _nodes[2];
		return new Tri(tri_nodes, _value);
	}

	// this should not happen
	return NULL;
}

} // end namespace MeshLib

