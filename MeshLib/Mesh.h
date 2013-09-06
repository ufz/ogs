/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESH_H_
#define MESH_H_

#include <cstdlib>
#include <string>
#include <vector>

#include "MeshEnums.h"

#include "BaseLib/Counter.h"

namespace MeshLib
{
	class Node;
	class Element;

/**
 * A basic mesh.
 */
class Mesh : BaseLib::Counter<Mesh>
{
	/* friend functions: */
	friend void removeMeshNodes(MeshLib::Mesh &mesh, const std::vector<std::size_t> &nodes);

public:
	/// Constructor using a mesh name and an array of nodes and elements
	Mesh(const std::string &name,
	     const std::vector<Node*> &nodes,
	     const std::vector<Element*> &elements);

	/// Copy constructor
	Mesh(const Mesh &mesh);

	/// Destructor
	virtual ~Mesh();

	/// Add a node to the mesh.
	void addNode(Node* node);

	/// Add an element to the mesh.
	void addElement(Element* elem);

	/// Returns the dimension of the mesh (determined by the maximum dimension over all elements).
	unsigned getDimension() const { return _mesh_dimension; }

	/// Get the node with the given index.
	const Node* getNode(unsigned idx) const { return _nodes[idx]; }

	/// Get the element with the given index.
	const Element* getElement(unsigned idx) const { return _elements[idx]; }

	/// Get the minimum edge length over all elements of the mesh.
	double getMinEdgeLength() const { return _edge_length[0]; }

	/// Get the maximum edge length over all elements of the mesh.
	double getMaxEdgeLength() const { return _edge_length[1]; }

	/// Get the number of elements
	std::size_t getNElements() const { return _elements.size(); }

	/// Get the number of nodes
	std::size_t getNNodes() const { return _nodes.size(); }

	/// Get name of the mesh.
	const std::string getName() const { return _name; }

	/// Get the nodes-vector for the mesh.
	std::vector<Node*> const& getNodes() const { return _nodes; }

	/// Get the element-vector for the mesh.
	std::vector<Element*> const& getElements() const { return _elements; }

	/// Resets the IDs of all mesh-elements to their position in the element vector
	void resetElementIDs();
	
	/// Resets the IDs of all mesh-nodes to their position in the node vector
	void resetNodeIDs();

	/// Changes the name of the mesh.
	void setName(const std::string &name) { this->_name = name; }

	/// Get id of the mesh
	std::size_t getID() const {return _id; }

protected:
	/// Set the minimum and maximum length over the edges of the mesh.
	void calcEdgeLengthRange();

	/**
	 * Resets the connected elements for the node vector, i.e. removes the old information and
	 * calls setElementsConnectedToNodes to set the new information.
	 * \attention This needs to be called if node neighbourhoods are reset.
	 */
	void resetElementsConnectedToNodes();	

	/// Sets the dimension of the mesh.
	void setDimension();

	/// Fills in the neighbor-information for nodes (i.e. which element each node belongs to).
	void setElementsConnectedToNodes();

	/// Fills in the neighbor-information for elements.
	void setNeighborElements();

	void setNodesConnectedByEdges();

	void setNodesConnectedByElements();

	std::size_t const _id;
	unsigned _mesh_dimension;
	double _edge_length[2];
	std::string _name;
	std::vector<Node*> _nodes;
	std::vector<Element*> _elements;

}; /* class */

} /* namespace */

#endif /* MESH_H_ */

