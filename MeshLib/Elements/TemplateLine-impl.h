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
                                                unsigned value, std::size_t id)
	: Element(value, id)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[2];
	std::fill_n(_neighbors, 2, nullptr);
	this->_length = this->computeVolume();
}

template<unsigned NNODES, CellType CELLLINETYPE>
TemplateLine<NNODES,CELLLINETYPE>::TemplateLine(Node* nodes[NNODES], unsigned value, std::size_t id)
	: Element(value, id)
{
	_nodes = nodes;
	_neighbors = new Element*[2];
	std::fill_n(_neighbors, 2, nullptr);
	this->_length = this->computeVolume();
}

template <unsigned NNODES, CellType CELLLINETYPE>
TemplateLine<NNODES,CELLLINETYPE>::TemplateLine(const TemplateLine<NNODES,CELLLINETYPE> &line)
	: Element(line.getValue(), line.getID())
{
	_nodes = new Node*[NNODES];
	for (unsigned k(0); k<NNODES; k++)
		_nodes[k] = line._nodes[k];

	_neighbors = new Element*[2];
	_neighbors[0] = line._neighbors[0];
	_neighbors[1] = line._neighbors[1];
	_length = line.getLength();
}

template <unsigned NNODES, CellType CELLLINETYPE>
TemplateLine<NNODES,CELLLINETYPE>::~TemplateLine()
{}

template <unsigned NNODES, CellType CELLLINETYPE>
ElementErrorCode TemplateLine<NNODES,CELLLINETYPE>::validate() const
{ 
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume] = this->hasZeroVolume();
	return error_code;
}

} // namespace MeshLib

