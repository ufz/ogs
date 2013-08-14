/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief Implementation of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Mesh.h"

#include "Node.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"

#include "logog/include/logog.hpp"

#include "RunTime.h"
#include "uniqueInsert.h"

namespace MeshLib
{

Mesh::Mesh(const std::string &name,
           const std::vector<Node*> &nodes,
           const std::vector<Element*> &elements)
	: _id(_counter_value), _mesh_dimension(0), _name(name), _nodes(nodes), _elements(elements)
{
	this->resetNodeIDs();
	this->resetElementIDs();
	this->setDimension();
	this->setElementsConnectedToNodes();
	//this->setNodesConnectedByEdges();
	//this->setNodesConnectedByElements();
	this->setElementNeighbors();

	_edge_length[0] =  std::numeric_limits<double>::max();
	_edge_length[1] = -std::numeric_limits<double>::max();
	this->calcEdgeLengthRange();
}

Mesh::Mesh(const Mesh &mesh)
	: _id(_counter_value), _mesh_dimension(mesh.getDimension()),
	  _name(mesh.getName()), _nodes(mesh.getNNodes()), _elements(mesh.getNElements())
{
	const std::vector<Node*> nodes (mesh.getNodes());
	const size_t nNodes (nodes.size());
	for (unsigned i=0; i<nNodes; ++i)
		_nodes[i] = new Node(*nodes[i]);

	const std::vector<Element*> elements (mesh.getElements());
	const size_t nElements (elements.size());
	for (unsigned i=0; i<nElements; ++i)
	{
		const size_t nElemNodes = elements[i]->getNNodes();
		_elements[i] = elements[i]->clone();
		for (unsigned j=0; j<nElemNodes; ++j)
			_elements[i]->_nodes[j] = _nodes[elements[i]->getNode(j)->getID()];
	}

	if (_mesh_dimension==0) this->setDimension();
	this->_edge_length[0] = mesh.getMinEdgeLength();
	this->_edge_length[1] = mesh.getMaxEdgeLength();
	this->setElementsConnectedToNodes();
	//this->setNodesConnectedByEdges();
	//this->setNodesConnectedByElements();
	this->setElementNeighbors();
}

Mesh::~Mesh()
{
	const size_t nElements (_elements.size());
	for (size_t i=0; i<nElements; ++i)
		delete _elements[i];

	const size_t nNodes (_nodes.size());
	for (size_t i=0; i<nNodes; ++i)
		delete _nodes[i];
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
	for (unsigned i=0; i<nNodes; ++i)
		elem->_nodes[i]->addElement(elem);
}

void Mesh::resetNodeIDs()
{
	const size_t nNodes (this->_nodes.size());
	for (unsigned i=0; i<nNodes; ++i)
		_nodes[i]->setID(i);
}

void Mesh::resetElementIDs()
{
	const size_t nElements (this->_elements.size());
	for (unsigned i=0; i<nElements; ++i)
		_elements[i]->setID(i);
}

void Mesh::setDimension()
{
	const size_t nElements (_elements.size());
	for (unsigned i=0; i<nElements; ++i)
		if (_elements[i]->getDimension() > _mesh_dimension)
			_mesh_dimension = _elements[i]->getDimension();
}

void Mesh::setElementsConnectedToNodes()
{
	const size_t nElements (_elements.size());
	for (unsigned i=0; i<nElements; ++i)
	{
		MeshLib::Element* element = _elements[i];
		const unsigned nNodes (element->getNNodes());
		for (unsigned j=0; j<nNodes; ++j)
			element->_nodes[j]->addElement(element);
	}
}

void Mesh::resetElementsConnectedToNodes()
{
	for (auto node = _nodes.begin(); node != _nodes.end(); ++node)
		if (*node)
			(*node)->_elements.clear();
	this->setElementsConnectedToNodes();
}

void Mesh::calcEdgeLengthRange()
{
	double min_length(0);
	double max_length(0);
	const std::size_t nElems (this->getNElements());
	for (std::size_t i=0; i<nElems; ++i)
	{
		_elements[i]->computeSqrEdgeLengthRange(min_length, max_length);
		this->_edge_length[0] = std::min(this->_edge_length[0], min_length);
		this->_edge_length[1] = std::max(this->_edge_length[1], max_length);
	}
	this->_edge_length[0] = sqrt(this->_edge_length[0]);
	this->_edge_length[1] = sqrt(this->_edge_length[1]);
}

void Mesh::setElementNeighbors()
{
	const size_t nElements = _elements.size();
	std::vector<Element*> neighbors;
	for (unsigned m(0); m<nElements; ++m)
	{
		// create vector with all elements connected to current element (includes lots of doubles!)
		Element *const element (_elements[m]);

		const size_t nNodes (element->getNNodes());
		for (unsigned n(0); n<nNodes; ++n)
		{
			std::vector<Element*> const& conn_elems ((element->getNode(n)->getElements()));
			neighbors.insert(neighbors.end(), conn_elems.begin(), conn_elems.end());
		}
			std::sort(neighbors.begin(), neighbors.end());
			auto const neighbors_new_end = std::unique(neighbors.begin(), neighbors.end());

			const unsigned nNeighbors ( std::distance(neighbors.begin(), neighbors_new_end));

		for (unsigned i(0); i<nNeighbors; ++i)
		{
			unsigned const opposite_face_id = element->addNeighbor(neighbors[i]);
			if (element->addNeighbor(neighbors[i]))
			{
				neighbors[i]->addNeighbor(element);
			}
		}
		neighbors.clear();
	}
}

void Mesh::setNodesConnectedByEdges()
{
	const size_t nNodes (this->_nodes.size());
	for (unsigned i=0; i<nNodes; ++i)
	{
		MeshLib::Node* node (_nodes[i]);
		std::vector<MeshLib::Node*> conn_set;
		const std::vector<MeshLib::Element*> &conn_elems (node->getElements());
		const size_t nConnElems (conn_elems.size());
		for (unsigned j=0; j<nConnElems; ++j)
		{
			const unsigned idx (conn_elems[j]->getNodeIDinElement(node));
			const unsigned nElemNodes (conn_elems[j]->getNNodes());
			for (unsigned k(0); k<nElemNodes; ++k)
			{
				bool is_in_vector (false);
				const size_t nConnNodes (conn_set.size());
				for (unsigned l(0); l<nConnNodes; ++l)
					if (conn_elems[j]->getNode(k) == conn_set[l])
						is_in_vector = true;
				if (is_in_vector) continue;
				if (conn_elems[j]->isEdge(idx, k))
					conn_set.push_back(_nodes[conn_elems[j]->getNode(k)->getID()]);
			}
		}
		node->setConnectedNodes(conn_set);
	}
}

void Mesh::setNodesConnectedByElements()
{
	const size_t nNodes (this->_nodes.size());
	for (unsigned i=0; i<nNodes; ++i)
	{
		MeshLib::Node* node (_nodes[i]);
		std::vector<MeshLib::Node*> conn_vec;
		const std::vector<MeshLib::Element*> &conn_elems (node->getElements());
		const size_t nConnElems (conn_elems.size());
		for (unsigned j=0; j<nConnElems; ++j)
		{
			const unsigned nElemNodes (conn_elems[j]->getNNodes());
			for (unsigned k(0); k<nElemNodes; ++k)
			{
				bool is_in_vector (false);
				const MeshLib::Node* c_node (conn_elems[j]->getNode(k));
				if (c_node == node) continue;
				const size_t nConnNodes (conn_vec.size());
				for (unsigned l(0); l<nConnNodes; ++l)
					if (c_node == conn_vec[l])
						is_in_vector = true;
				if (!is_in_vector)
					conn_vec.push_back(_nodes[c_node->getID()]);
			}
		}
		node->setConnectedNodes(conn_vec);
	}
}


}
