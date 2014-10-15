/**
 * \file   ElementStatus.cpp
 * \author Karsten Rink
 * \date   2012-12-18
 * \brief  Implementation of the ElementStatus class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementStatus.h"

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

namespace MeshLib {

ElementStatus::ElementStatus(Mesh const*const mesh)
: _mesh(mesh), _element_status(mesh->getNElements(), true)
{
	const std::vector<MeshLib::Node*> &nodes (_mesh->getNodes());
	for (auto node = nodes.cbegin(); node != nodes.cend(); ++node)
		_active_nodes.push_back((*node)->getNElements());
}


std::vector<std::size_t> ElementStatus::getActiveElements() const
{
	std::vector<std::size_t> active_elements;
	active_elements.reserve(this->getNActiveElements());
	const std::size_t nElems (_mesh->getNElements());
	for (std::size_t i=0; i<nElems; ++i)
		if (_element_status[i])
			active_elements.push_back(i);
	return active_elements;
}

std::vector<std::size_t> ElementStatus::getActiveNodes() const 
{
	std::vector<std::size_t> active_nodes;
	active_nodes.reserve(this->getNActiveNodes());
	const std::size_t nNodes (_mesh->getNNodes());
	for (std::size_t i=0; i<nNodes; ++i)
		if (_active_nodes[i]>0)
			active_nodes.push_back(i);
	return active_nodes;
}

std::vector<std::size_t> ElementStatus::getActiveElementsAtNode(std::size_t node_id) const
{
	const std::size_t nActiveElements (_active_nodes[node_id]);
	const std::vector<Element*> &elements (_mesh->getNode(node_id)->getElements());
	std::vector<std::size_t> active_elements;
	active_elements.reserve(nActiveElements);
	for (auto elem = elements.cbegin(); elem != elements.cend(); ++elem)
	{
		if (active_elements.size() == nActiveElements)
			return active_elements;
		const std::size_t id ((*elem)->getID());
		if (_element_status[id])
			active_elements.push_back(id);
	}
	return active_elements;
}

std::size_t ElementStatus::getNActiveNodes() const 
{
	return _active_nodes.size() - std::count(_active_nodes.cbegin(), _active_nodes.cend(), 0);
}

std::size_t ElementStatus::getNActiveElements() const 
{
	return static_cast<std::size_t>(std::count(_element_status.cbegin(), _element_status.cend(), true));
}

void ElementStatus::setAll(bool status)
{
	std::fill(_element_status.begin(), _element_status.end(), status);

	if (status)
	{
		const std::vector<MeshLib::Node*> &nodes (_mesh->getNodes());
		const std::size_t nNodes (_mesh->getNNodes());
		for (std::size_t i=0; i<nNodes; ++i)
			_active_nodes[i] = nodes[i]->getNElements();
	}
	else
		std::fill(_active_nodes.begin(), _active_nodes.end(), 0);
}

void ElementStatus::setElementStatus(std::size_t i, bool status)
{
	if (_element_status[i] != status)
	{
		const int change = (status) ? 1 : -1;
		_element_status[i] = status;
		const unsigned nElemNodes (_mesh->getElement(i)->getNBaseNodes());
		MeshLib::Node const*const*const nodes = _mesh->getElement(i)->getNodes();
		for (unsigned i=0; i<nElemNodes; ++i)
		{
			assert(_active_nodes[i]<255); // if one node has >255 connected elements the data type is too small
			_active_nodes[nodes[i]->getID()] += change;
		}
	}
}

void ElementStatus::setMaterialStatus(unsigned material_id, bool status)
{
	const std::size_t nElems (_mesh->getNElements());
	for (std::size_t i=0; i<nElems; ++i)
		if (_mesh->getElement(i)->getValue() == material_id)
			this->setElementStatus(i, status);
}

}

