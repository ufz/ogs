/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad2<NNODES, CELLQUADTYPE>::n_all_nodes;

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad2<NNODES, CELLQUADTYPE>::n_base_nodes;

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad2<NNODES, CELLQUADTYPE>::_edge_nodes[4][3] =
{
	{0, 1, 4}, // Edge 0
	{1, 2, 5}, // Edge 1
	{2, 3, 6}, // Edge 2
	{0, 3, 7}  // Edge 3
};

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad2<NNODES,CELLQUADTYPE>::TemplateQuad2(Node* nodes[NNODES], unsigned value, std::size_t id)
	: Face(value, id)
{
	_nodes = nodes;

	_neighbors = new Element*[4];
	std::fill(_neighbors, _neighbors + 4, nullptr);

	this->_area = this->computeVolume();
}

template<unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad2<NNODES,CELLQUADTYPE>::TemplateQuad2(std::array<Node*, NNODES> const& nodes,
                                                unsigned value, std::size_t id)
	: Face(value, id)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[4];
	std::fill(_neighbors, _neighbors + 4, nullptr);

	this->_area = this->computeVolume();
}

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad2<NNODES,CELLQUADTYPE>::TemplateQuad2(const TemplateQuad2<NNODES,CELLQUADTYPE> &quad)
	: Face(quad.getValue(), quad.getID())
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

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad2<NNODES,CELLQUADTYPE>::~TemplateQuad2()
{
}

template <unsigned NNODES, CellType CELLQUADTYPE>
double TemplateQuad2<NNODES,CELLQUADTYPE>::computeVolume()
{
	//TODO
	return GeoLib::calcTriangleArea(*_nodes[0], *_nodes[1], *_nodes[2])
         + GeoLib::calcTriangleArea(*_nodes[2], *_nodes[3], *_nodes[0]);
}

template <unsigned NNODES, CellType CELLQUADTYPE>
bool TemplateQuad2<NNODES,CELLQUADTYPE>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<4; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType CELLQUADTYPE>
bool TemplateQuad2<NNODES,CELLQUADTYPE>::isPntInElement(GeoLib::Point const& pnt, double eps) const
{
	//TODO
    bool a (GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[1], *_nodes[2], eps));
    bool b ( GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[2], *_nodes[3], eps));
	return (GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[1], *_nodes[2], eps) ||
	        GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[2], *_nodes[3], eps));
}

template <unsigned NNODES, CellType CELLQUADTYPE>
Element* TemplateQuad2<NNODES,CELLQUADTYPE>::clone() const
{
	return new TemplateQuad2(*this);
}

template <unsigned NNODES, CellType CELLQUADTYPE>
unsigned TemplateQuad2<NNODES,CELLQUADTYPE>::identifyFace(Node* nodes[3]) const
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

template <unsigned NNODES, CellType CELLQUADTYPE>
ElementErrorCode TemplateQuad2<NNODES,CELLQUADTYPE>::validate() const
{
	//TODO
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume]  = this->hasZeroVolume();
	error_code[ElementErrorFlag::NonCoplanar] = (!GeoLib::isCoplanar(*_nodes[0], *_nodes[1], *_nodes[2], *_nodes[3]));
	// for collapsed quads (i.e. reduced to a line) this test might result "false" as all four points are actually located on a line.
	if (!error_code[ElementErrorFlag::ZeroVolume]) 
		error_code[ElementErrorFlag::NonConvex]   = (!(GeoLib::dividedByPlane(*_nodes[0], *_nodes[2], *_nodes[1], *_nodes[3]) &&
			                                           GeoLib::dividedByPlane(*_nodes[1], *_nodes[3], *_nodes[0], *_nodes[2])));
	error_code[ElementErrorFlag::NodeOrder]  = !this->testElementNodeOrder();
	return error_code;
}

}

