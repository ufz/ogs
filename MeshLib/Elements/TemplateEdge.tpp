/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Implementation of the TemplateEdge class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace MeshLib {

template<unsigned NNODES, CellType::type CELLEDGETYPE>
TemplateEdge<NNODES,CELLEDGETYPE>::TemplateEdge(Node* nodes[NNODES], unsigned value) :
	Element(value)
{
	_nodes = nodes;
	this->_length = this->computeVolume();
}

template <unsigned NNODES, CellType::type CELLEDGETYPE>
TemplateEdge<NNODES,CELLEDGETYPE>::TemplateEdge(const TemplateEdge<NNODES,CELLEDGETYPE> &edge) :
	Element(edge.getValue())
{
	_nodes = new Node*[NNODES];
	for (unsigned k(0); k<NNODES; k++)
		_nodes[k] = edge._nodes[k];
	_length = edge.getLength();
}

template <unsigned NNODES, CellType::type CELLEDGETYPE>
TemplateEdge<NNODES,CELLEDGETYPE>::~TemplateEdge()
{}

} // namespace MeshLib

