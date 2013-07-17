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
	this->resetNodeIDs(); // reset node ids so they match the node position in the vector
	_edge_length[0] = 0;
	_edge_length[1] = 0;
	this->setDimension();
	this->makeNodesUnique();
	this->setElementsConnectedToNodes();
	//this->setNodesConnectedByEdges();
	//this->setNodesConnectedByElements();
	this->setElementsConnectedToElements();
	this->removeUnusedMeshNodes();
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
	this->setElementsConnectedToNodes();
	//this->setNodesConnectedByEdges();
	//this->setNodesConnectedByElements();
	this->setElementsConnectedToElements();
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

void Mesh::makeNodesUnique()
{
	//check for unique mesh nodes
	//PointVec::makePntsUnique

	//replace node pointers in elements
	unsigned nElements (_elements.size());
	for (unsigned i=0; i<nElements; ++i)
	{
		unsigned nNodes (_elements[i]->getNNodes());
		for (unsigned j=0; j<nNodes; ++j)
			_elements[i]->getNodeIndex(j);
	}

	//set correct id for each node

	//if (this->getDimension() > 1)
	//	this->removeMeshElements(MshElemType::EDGE);

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
//#ifndef NDEBUG
	// search for nodes that are not part of any element
	unsigned count(0);
	const size_t nNodes (_nodes.size());
	for (unsigned i=0; i<nNodes; ++i)
		if (_nodes[i]->getNElements() == 0)
		{
			WARN ("Node %d is not part of any element.", i);
			++count;
		}
	if (count)
		WARN ("%d unused mesh nodes found.", count);
//#endif
}

void Mesh::setEdgeLengthRange(const double &min_length, const double &max_length)
{
	if (min_length <= max_length)
	{
		_edge_length[0] = min_length;
		_edge_length[1] = max_length;
	}
	else
		ERR ("Error in MeshLib::Mesh::setEdgeLengthRange() - min length > max length.");
}

void Mesh::setElementsConnectedToElements()
{
	const size_t nElements = _elements.size();
	for (unsigned m(0); m<nElements; ++m)
	{
		// create vector with all elements connected to current element (includes lots of doubles!)
		std::vector<Element*> neighbors;
		Element *const element (_elements[m]);
		if (element->getGeomType() != MshElemType::EDGE)
		{
			const size_t nNodes (element->getNNodes());
			for (unsigned n(0); n<nNodes; ++n)
			{
				std::vector<Element*> const& conn_elems ((element->getNode(n)->getElements()));
				neighbors.insert(neighbors.end(), conn_elems.begin(), conn_elems.end());
			}

			const unsigned nNeighbors ( neighbors.size() );

			for (unsigned i(0); i<nNeighbors; ++i)
			{
				if (element->addNeighbor(neighbors[i]) && neighbors[i]->getGeomType() != MshElemType::EDGE)
				{
					neighbors[i]->addNeighbor(element);
				}
			}
		}
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

void Mesh::addDoublePropertyVec(std::string const& prop_name, std::vector<double> const& elem_props)
{
	_double_prop_vecs.push_back(std::pair<std::string,std::vector<double> >(prop_name, elem_props));
}

void Mesh::addUnsignedPropertyVec(std::string const& prop_name, std::vector<unsigned> const& elem_props)
{
	_unsigned_prop_vecs.push_back(std::pair<std::string,std::vector<unsigned> >(prop_name, elem_props));
}

boost::optional<std::vector<double> const&>
Mesh::getDoublePropertyVec(std::string const& prop_name) const
{
	auto it = _double_prop_vecs.begin();
	while (it != _double_prop_vecs.end() && it->first.compare(prop_name) != 0)
		++it;

	if (it == _double_prop_vecs.end())
		return boost::optional<std::vector<double> const&>();

	return boost::optional<std::vector<double> const&>(it->second);
}

boost::optional<std::vector<unsigned> const&>
Mesh::getUnsignedPropertyVec(std::string const& prop_name) const
{
	auto it = _unsigned_prop_vecs.begin();
	while (it != _unsigned_prop_vecs.end() && it->first.compare(prop_name) != 0)
		++it;

	if (it == _unsigned_prop_vecs.end())
		return boost::optional<std::vector<unsigned> const&>();

	return boost::optional<std::vector<unsigned> const&>(it->second);
}

void Mesh::removeUnusedMeshNodes()
{
	unsigned count(0);
	std::vector<MeshLib::Node*>::iterator it (this->_nodes.begin());
	while(it != this->_nodes.end())
	{
		if ((*it)->getNElements() == 0)
		{
			delete *it;
			*it = nullptr;
			++it;
			++count;
		}
		else ++it;
	}
	auto node_vec_end = std::remove(_nodes.begin(), _nodes.end(), nullptr);
	_nodes.erase(node_vec_end, _nodes.end());

	if (count)
	{
		INFO("Removed %d unused mesh nodes.", count );
		this->resetNodeIDs();
	}
}

void Mesh::removeMeshElements(MshElemType::type t)
{
	unsigned count(0);
	for (std::vector<MeshLib::Element*>::iterator it = this->_elements.begin(); it != this->_elements.end();)
	{
		if ((*it)->getGeomType() == t)
		{
			delete *it;
			it = this->_elements.erase(it);
			++count;
		}
		else
			++it;
	}
	INFO("Removed %d  elements of type %s from mesh.", count, MshElemType2String(t).c_str());
}
}
