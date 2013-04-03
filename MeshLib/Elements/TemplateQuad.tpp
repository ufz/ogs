/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Quad class-
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <array>

#include "Node.h"
#include "Tri.h"

// MathLib
#include "MathTools.h"
#include "AnalyticalGeometry.h"

namespace MeshLib
{

template <unsigned NNODES, CellType::type CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::TemplateQuad(Node* nodes[NNODES], unsigned value)
	: Face(value)
{
	_nodes = nodes;

	_neighbors = new Element*[4];
	std::fill(_neighbors, _neighbors + 4, nullptr);

	this->_area = this->computeVolume();
}

template<unsigned NNODES, CellType::type CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::TemplateQuad(std::array<Node*, NNODES> const& nodes,
                                                unsigned value)
	: Face(value)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[4];
	std::fill(_neighbors, _neighbors + 4, nullptr);

	this->_area = this->computeVolume();
}

template <unsigned NNODES, CellType::type CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::TemplateQuad(const TemplateQuad<NNODES,CELLQUADTYPE> &quad)
	: Face(quad.getValue())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++) {
		_nodes[i] = quad._nodes[i];
	}

	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++) {
		_neighbors[i] = quad._neighbors[i];
	}

	_area = quad.getArea();
}

template <unsigned NNODES, CellType::type CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::~TemplateQuad()
{
}

template <unsigned NNODES, CellType::type CELLQUADTYPE>
double TemplateQuad<NNODES,CELLQUADTYPE>::computeVolume()
{
	return MathLib::calcTriangleArea(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords())
         + MathLib::calcTriangleArea(_nodes[2]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords());
}

template <unsigned NNODES, CellType::type CELLQUADTYPE>
bool TemplateQuad<NNODES,CELLQUADTYPE>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<4; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType::type CELLQUADTYPE>
bool TemplateQuad<NNODES,CELLQUADTYPE>::isPntInside(GeoLib::Point const& pnt, double eps) const
{
	return (GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[1], *_nodes[2], eps) ||
					GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[2], *_nodes[3], eps));
}

template <unsigned NNODES, CellType::type CELLQUADTYPE>
Element* TemplateQuad<NNODES,CELLQUADTYPE>::clone() const
{
	return new TemplateQuad(*this);
}

template <unsigned NNODES, CellType::type CELLQUADTYPE>
unsigned TemplateQuad<NNODES,CELLQUADTYPE>::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<4; i++)
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

template <unsigned NNODES, CellType::type CELLQUADTYPE>
Element* TemplateQuad<NNODES,CELLQUADTYPE>::reviseElement() const
{
	if (_nodes[0] == _nodes[1] || _nodes[1] == _nodes[2]) {
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		tri_nodes[0] = _nodes[0];
		tri_nodes[1] = _nodes[2];
		tri_nodes[2] = _nodes[3];
		return new Tri(tri_nodes, _value);
	}

	if (_nodes[2] == _nodes[3] || _nodes[3] == _nodes[0]) {
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		tri_nodes[0] = _nodes[0];
		tri_nodes[1] = _nodes[1];
		tri_nodes[2] = _nodes[2];
		return new Tri(tri_nodes, _value);
	}

	// this should not happen
	return NULL;
}

}

