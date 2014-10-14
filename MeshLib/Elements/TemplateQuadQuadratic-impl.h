/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include "MeshLib/Node.h"
#include "Line.h"

namespace MeshLib
{

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuadQuadratic<NNODES,CELLQUADTYPE>::TemplateQuadQuadratic(Node* nodes[NNODES], unsigned value, std::size_t id)
	: TemplateQuad<NNODES, CELLQUADTYPE, detail::QuadEdgeQuadraticNodes>(nodes, value, id)
{}

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuadQuadratic<NNODES,CELLQUADTYPE>::TemplateQuadQuadratic(std::array<Node*, NNODES> const& nodes,
                                                unsigned value, std::size_t id)
	: TemplateQuad<NNODES, CELLQUADTYPE, detail::QuadEdgeQuadraticNodes>(nodes, value, id)
{}

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuadQuadratic<NNODES,CELLQUADTYPE>::TemplateQuadQuadratic(const TemplateQuadQuadratic<NNODES,CELLQUADTYPE> &quad)
	: TemplateQuad<NNODES, CELLQUADTYPE, detail::QuadEdgeQuadraticNodes>(quad)
{}

template <unsigned NNODES, CellType CELLQUADTYPE>
const Element* TemplateQuadQuadratic<NNODES,CELLQUADTYPE>::getEdge(unsigned i) const
{
	if (i < this->getNEdges())
	{
		Node** nodes = new Node*[3];
		nodes[0] = const_cast<Node*>(this->getEdgeNode(i,0));
		nodes[1] = const_cast<Node*>(this->getEdgeNode(i,1));
		nodes[2] = const_cast<Node*>(this->getEdgeNode(i,2));
		return new Line(nodes);
	}
	ERR("Error in MeshLib::Element::getEdge() - Index does not exist.");
	return nullptr;
}

template <unsigned NNODES, CellType CELLQUADTYPE>
Element* TemplateQuadQuadratic<NNODES,CELLQUADTYPE>::clone() const
{
	return new TemplateQuadQuadratic(*this);
}

}

