/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the TemplateTet class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include "Node.h"
#include "Tri.h"

#include "MathTools.h"

namespace MeshLib {

template <unsigned NNODES, CellType::type CELLTETTYPE>
const unsigned TemplateTet<NNODES,CELLTETTYPE>::_face_nodes[4][3] =
{
	{0, 2, 1}, // Face 0
	{0, 1, 3}, // Face 1
	{1, 2, 3}, // Face 2
	{2, 0, 3}  // Face 3
};

template <unsigned NNODES, CellType::type CELLTETTYPE>
const unsigned TemplateTet<NNODES,CELLTETTYPE>::_edge_nodes[6][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{0, 2}, // Edge 2
	{0, 3}, // Edge 3
	{1, 3}, // Edge 4
	{2, 3}  // Edge 5
};

template <unsigned NNODES, CellType::type CELLTETTYPE>
TemplateTet<NNODES,CELLTETTYPE>::TemplateTet(Node* nodes[NNODES], unsigned value)
	: Cell(value)
{
	_nodes = nodes;

	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++)
		_neighbors[i] = NULL;

	this->_volume = this->computeVolume();
}

template<unsigned NNODES, CellType::type CELLTETTYPE>
TemplateTet<NNODES,CELLTETTYPE>::TemplateTet(std::array<Node*, NNODES> const& nodes,
                                             unsigned value)
	: Cell(value)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[4];
	std::fill(_neighbors, _neighbors + 4, nullptr);

	this->_volume = this->computeVolume();
}

template <unsigned NNODES, CellType::type CELLTETTYPE>
TemplateTet<NNODES,CELLTETTYPE>::TemplateTet(const TemplateTet<NNODES,CELLTETTYPE> &tet)
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

template <unsigned NNODES, CellType::type CELLTETTYPE>
TemplateTet<NNODES,CELLTETTYPE>::~TemplateTet()
{
}

template <unsigned NNODES, CellType::type CELLTETTYPE>
double TemplateTet<NNODES,CELLTETTYPE>::computeVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[3]->getCoords());
}

template <unsigned NNODES, CellType::type CELLTETTYPE>
const Element* TemplateTet<NNODES,CELLTETTYPE>::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];
		return new Tri(nodes);
	}
	ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
	return NULL;
}

template <unsigned NNODES, CellType::type CELLTETTYPE>
bool TemplateTet<NNODES,CELLTETTYPE>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<6; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType::type CELLTETTYPE>
Element* TemplateTet<NNODES,CELLTETTYPE>::clone() const
{
	return new TemplateTet<NNODES,CELLTETTYPE>(*this);
}

template <unsigned NNODES, CellType::type CELLTETTYPE>
unsigned TemplateTet<NNODES,CELLTETTYPE>::identifyFace(Node* nodes[3]) const
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

template <unsigned NNODES, CellType::type CELLTETTYPE>
Element* TemplateTet<NNODES,CELLTETTYPE>::reviseElement() const
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

