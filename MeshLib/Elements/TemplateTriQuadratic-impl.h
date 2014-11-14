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
TemplateTriQuadratic<NNODES,CELLQUADTYPE>::TemplateTriQuadratic(Node* nodes[NNODES], unsigned value, std::size_t id)
	: TemplateTri<NNODES, CELLQUADTYPE, detail::TriEdgeQuadraticNodes>(nodes, value, id)
{}

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateTriQuadratic<NNODES,CELLQUADTYPE>::TemplateTriQuadratic(std::array<Node*, NNODES> const& nodes,
                                                unsigned value, std::size_t id)
	: TemplateTri<NNODES, CELLQUADTYPE, detail::TriEdgeQuadraticNodes>(nodes, value, id)
{}

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateTriQuadratic<NNODES,CELLQUADTYPE>::TemplateTriQuadratic(const TemplateTriQuadratic<NNODES,CELLQUADTYPE> &tri)
	: TemplateTri<NNODES, CELLQUADTYPE, detail::TriEdgeQuadraticNodes>(tri)
{}

template <unsigned NNODES, CellType CELLQUADTYPE>
const Element* TemplateTriQuadratic<NNODES,CELLQUADTYPE>::getEdge(unsigned i) const
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
Element* TemplateTriQuadratic<NNODES,CELLQUADTYPE>::clone() const
{
	return new TemplateTriQuadratic(*this);
}

}

