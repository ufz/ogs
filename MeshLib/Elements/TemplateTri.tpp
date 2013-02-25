/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Implementation of the TemplateTri class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// MathLib
#include "AnalyticalGeometry.h"

namespace MeshLib {

template <unsigned NNODES, CellType::type CELLTRITYPE>
TemplateTri<NNODES,CELLTRITYPE>::TemplateTri(Node* nodes[NNODES], unsigned value) :
	Face(value)
{
	_nodes = nodes;
	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++)
		_neighbors[i] = NULL;
	this->_area = this->computeVolume();
}

template <unsigned NNODES, CellType::type CELLTRITYPE>
TemplateTri<NNODES,CELLTRITYPE>::TemplateTri(const TemplateTri<NNODES,CELLTRITYPE> &tri) :
	Face(tri.getValue())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++)
	{
		_nodes[i] = tri._nodes[i];
	}

	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++) {
		_neighbors[i] = tri._neighbors[i];
	}

	_area = tri.getArea();
}

template <unsigned NNODES, CellType::type CELLTRITYPE>
TemplateTri<NNODES,CELLTRITYPE>::~TemplateTri()
{}

template <unsigned NNODES, CellType::type CELLTRITYPE>
bool TemplateTri<NNODES,CELLTRITYPE>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<3; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType::type CELLTRITYPE>
bool TemplateTri<NNODES,CELLTRITYPE>::isPntInside(GeoLib::Point const& pnt, double eps) const
{
	return GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[1], *_nodes[2], eps);
}

template <unsigned NNODES, CellType::type CELLTRITYPE>
Element* TemplateTri<NNODES,CELLTRITYPE>::reviseElement() const
{
	// try to create an edge
	if (_nodes[0] == _nodes[1] || _nodes[1] == _nodes[2]) {
		Node** nodes = new Node*[2];
		nodes[0] = _nodes[0];
		nodes[1] = _nodes[2];
		return new Edge(nodes, _value);
	}

	if (_nodes[0] == _nodes[2]) {
		Node** nodes = new Node*[2];
		nodes[0] = _nodes[0];
		nodes[1] = _nodes[1];
		return new Edge(nodes, _value);
	}

	return NULL;
}

template <unsigned NNODES, CellType::type CELLTRITYPE>
unsigned TemplateTri<NNODES,CELLTRITYPE>::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<3; i++)
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

} // namespace MeshLib

