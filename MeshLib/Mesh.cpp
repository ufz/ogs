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
	: _mesh_dimension(0), _name(name), _nodes(nodes), _elements(elements)
{
	this->resetNodeIDs(); // reset node ids so they match the node position in the vector
	_edge_length[0] = 0;
	_edge_length[1] = 0;
	this->makeNodesUnique();
	this->setDimension();
	this->setElementInformationForNodes();
	this->setNeighborInformationForElements();
}

Mesh::Mesh(const Mesh &mesh)
	: _mesh_dimension(mesh.getDimension()), _name(mesh.getName()), _nodes(mesh.getNodes()), _elements(mesh.getElements())
{
	const std::vector<Node*> nodes (mesh.getNodes());
	const size_t nNodes (nodes.size());
	for (unsigned i=0; i<nNodes; i++)
		_nodes[i] = new Node(*nodes[i]);

	const std::vector<Element*> elements (mesh.getElements());
	const size_t nElements (elements.size());
	for (unsigned i=0; i<nElements; i++)
	{
		const size_t nElemNodes = elements[i]->getNNodes();
		_elements[i] = elements[i]->clone();
		for (unsigned j=0; j<nElemNodes; j++)
			_elements[i]->_nodes[j] = _nodes[elements[i]->getNode(j)->getID()];
	}

	if (_mesh_dimension==0) this->setDimension();
	this->setElementInformationForNodes();
	this->setNeighborInformationForElements();
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

void Mesh::resetNodeIDs()
{
	const size_t nNodes (this->_nodes.size());
	for (unsigned i=0; i<nNodes; i++)
		_nodes[i]->setID(i);
}

void Mesh::setDimension()
{
	const size_t nElements (_elements.size());
	for (unsigned i=0; i<nElements; i++)
		if (_elements[i]->getDimension() > _mesh_dimension)
			_mesh_dimension = _elements[i]->getDimension();
}

void Mesh::setElementInformationForNodes()
{
	const size_t nElements (_elements.size());
#ifdef _OPENMP
	OPENMP_LOOP_TYPE i;
	#pragma omp parallel for
#else
	unsigned i(0)
#endif
	for (i=0; i<nElements; i++)
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
		std::cerr << "Error in MeshLib::Mesh::setEdgeLengthRange() - min length > max length." << std::endl;
}

void Mesh::setNeighborInformationForElements()
{
	const size_t nElements = _elements.size();
#ifdef _OPENMP
	OPENMP_LOOP_TYPE m;
	#pragma omp parallel for
#else
	unsigned m(0);
#endif
	for (m=0; m<nElements; m++)
	{
		// create vector with all elements connected to current element (includes lots of doubles!)
		std::vector<Element*> neighbors;
		Element *const element (_elements[m]);
		if (element->getType() != MshElemType::EDGE)
		{
			const size_t nNodes (element->getNNodes());
			for (unsigned n(0); n<nNodes; n++)
			{
				std::vector<Element*> const& conn_elems ((element->getNode(n)->getElements()));
				neighbors.insert(neighbors.end(), conn_elems.begin(), conn_elems.end());
			}

			const unsigned nNeighbors ( neighbors.size() );

			for (unsigned i(0); i<nNeighbors; i++)
			{
				if (element->addNeighbor(neighbors[i]) && neighbors[i]->getType() != MshElemType::EDGE)
				{
					neighbors[i]->addNeighbor(element);
				}
			}
		}
	}
}

}

