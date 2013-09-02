/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Implementation of the TemplateLine class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace MeshLib
{
template<unsigned NNODES, CellType CELLLINETYPE>
TemplateLine<NNODES,CELLLINETYPE>::TemplateLine(std::array<Node*, NNODES> const& nodes,
                                                unsigned value)
	: Edge(value)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	this->_length = this->computeVolume();
}

template<unsigned NNODES, CellType CELLLINETYPE>
TemplateLine<NNODES,CELLLINETYPE>::TemplateLine(Node* nodes[NNODES], unsigned value) :
	Edge(value)
{
	_nodes = nodes;
	this->_length = this->computeVolume();
}

template <unsigned NNODES, CellType CELLLINETYPE>
TemplateLine<NNODES,CELLLINETYPE>::TemplateLine(const TemplateLine<NNODES,CELLLINETYPE> &line) :
	Edge(line.getValue())
{
	_nodes = new Node*[NNODES];
	for (unsigned k(0); k<NNODES; k++)
		_nodes[k] = line._nodes[k];
	_length = line.getLength();
}

template <unsigned NNODES, CellType CELLLINETYPE>
TemplateLine<NNODES,CELLLINETYPE>::~TemplateLine()
{}

} // namespace MeshLib

