/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Mesh.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Mesh.h"

#include "Node.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"

namespace MeshLib {

Mesh::Mesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements)
	: _name(name), _nodes(nodes), _elements(elements)
{
	_edge_length[0] = 0;
	_edge_length[1] = 0;
	this->makeNodesUnique();
	this->setElementInformationForNodes();
	this->setNeighborInformationForElements();
}

Mesh::Mesh(const Mesh &mesh)
	: _name(mesh.getName()), _nodes(mesh.getNodes()), _elements(mesh.getElements())
{
}

Mesh::~Mesh()
{
	const size_t nElements (_elements.size());
	for (size_t i=0; i<nElements; i++)
		delete _elements[i];

	const size_t nNodes (_nodes.size());
	for (size_t i=0; i<nNodes; i++)
		delete _nodes[i];
}

void Mesh::makeNodesUnique()
{
	//check for unique mesh nodes
	//PointVec::makePntsUnique

	//replace node pointers in elements
	unsigned nElements (_elements.size());
	for (unsigned i=0; i<nElements; i++)
	{
		unsigned nNodes (_elements[i]->getNNodes());
		for (unsigned j=0; j<nNodes; j++)
			_elements[i]->getNodeIndex(j);
	}

	//set correct id for each node

}

void Mesh::addNode(Node* node)
{
	_nodes.push_back(node);
}

void Mesh::addElement(Element* elem)
{
	_elements.push_back(elem);

	// add element information to nodes
	unsigned nNodes (elem->getNNodes());
	for (unsigned i=0; i<nNodes; i++)
		elem->_nodes[i]->addElement(elem);
}

void Mesh::setElementInformationForNodes()
{
	const size_t nElements (_elements.size());
	for (unsigned i=0; i<nElements; i++)
	{
		const unsigned nNodes (_elements[i]->getNNodes());
		for (unsigned j=0; j<nNodes; j++)
			_elements[i]->_nodes[j]->addElement(_elements[i]);
	}
}

void Mesh::setEdgeLengthRange(const double &min_length, const double &max_length)
{
	if (min_length <= max_length)
	{
		_edge_length[0] = min_length;
		_edge_length[1] = max_length;
	}
	else
		std::cerr << "Error in MeshLib::Mesh::setEdgeLengthRange() - min length < max length." << std::endl;
}

void Mesh::setNeighborInformationForElements()
{
	/* TODO
	const size_t nElements(_elements.size());
	std::vector<std::vector<char>> nb (nElements, std::vector<char>(nElements));

	for (unsigned i=0; i<nElements; i++)
	{
		Element* elem = _elements[i];
		const size_t nNodes (elem->getNNodes());
		for (unsigned j=0; j<nNodes; j++)
		{
			const Node* node = elem->getNode(j);

		}
	}
	*/
}

}

