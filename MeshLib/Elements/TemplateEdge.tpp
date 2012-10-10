/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file TemplateEdge.tpp
 *
 *  Created on  Sep 27, 2012 by Thomas Fischer
 */

namespace MeshLib {

template <unsigned ORDER, unsigned NNODES>
TemplateEdge<ORDER,NNODES>::TemplateEdge(Node* nodes[NNODES], unsigned value) :
	Element(value)
{
	_nodes = nodes;
	this->_length = this->computeVolume();
}

template <unsigned ORDER, unsigned NNODES>
TemplateEdge<ORDER,NNODES>::TemplateEdge(const TemplateEdge<ORDER, NNODES> &edge) :
	Element(edge.getValue())
{
	_nodes = new Node*[NNODES];
	for (unsigned k(0); k<NNODES; k++)
		_nodes[k] = edge._nodes[k];
	_length = edge.getLength();
}

template <unsigned ORDER, unsigned NNODES>
TemplateEdge<ORDER,NNODES>::~TemplateEdge()
{}

} // namespace MeshLib

