/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Mesh.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef MESH_H_
#define MESH_H_

#include <cstdlib>
#include <string>
#include <vector>

#include "MshEnums.h"

namespace MeshLib
{
	class Node;
	class Element;

/**
 * A basic mesh.
 */
class Mesh
{

public:
	/// Constructor using a mesh name and an array of nodes and elements
	Mesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements);

	/// Copy constructor
	Mesh(const Mesh &mesh);

	/// Destructor
	virtual ~Mesh();

	/// Add a node to the mesh.
	void addNode(Node* node);

	/// Add an element to the mesh.
	void addElement(Element* elem);

	/// Returns the dimension of the mesh (determined by the maximum dimension over all elements).
	unsigned getDimension() const { return _mesh_dimension; };

	/// Get the minimum edge length over all elements of the mesh.
	double getMinEdgeLength() { return _edge_length[0]; };

	/// Get the maximum edge length over all elements of the mesh.
	double getMaxEdgeLength() { return _edge_length[1]; };

	/// Get the node with the given index.
	const Node* getNode(unsigned idx) const { return _nodes[idx]; };

	/// Get the element with the given index.
	const Element* getElement(unsigned idx) const { return _elements[idx]; };

	/// Get the minimum edge length for the mesh
	double getMinEdgeLength() const { return _edge_length[0]; };

	/// Get the maximum edge length for the mesh
	double getMaxEdgeLength() const { return _edge_length[1]; };

	/// Get the number of elements
	std::size_t getNElements() const { return _elements.size(); };

	/// Get the number of nodes
	std::size_t getNNodes() const { return _nodes.size(); };

	/// Get name of the mesh.
	const std::string getName() const { return _name; };

	/// Get the nodes-vector for the mesh.
	std::vector<Node*> const& getNodes() const { return _nodes; };

	/// Get the element-vector for the mesh.
	std::vector<Element*> const& getElements() const { return _elements; };

	/// Resets the IDs of all mesh-nodes to their position in the node vector
	void resetNodeIDs();

	/**
	 * Set the minimum and maximum length over the edges of the mesh.
	 * This should have been previously calcumlated using the Element::computeSqrEdgeLengthRange(min, max)
	 * function or by some other means.
	 */
	void setEdgeLengthRange(const double &min_length, const double &max_length);

	/// Changes the name of the mesh.
	void setName(const std::string &name) { this->_name = name; };

protected:
	/// Checks the coordinates of all mesh nodes and removes identical nodes. Elements are adapted accordingly.
	void makeNodesUnique();

	/// Removes nodes that are not part of any element.
	void removeUnusedMeshNodes();

	/// Removes elements of the given type t from a mesh
	void removeMeshElements(MshElemType::type t);

	/// Sets the dimension of the mesh.
	void setDimension();

	/// Fills in the neighbor-information for nodes (i.e. which element each node belongs to).
	void setElementsConnectedToNodes();

	/// Fills in the neighbor-information for elements.
	void setElementsConnectedToElements();

	void setNodesConnectedByEdges();

	void setNodesConnectedByElements();

	static std::vector<Node*> deepCopyNodes(Mesh* const*const mesh);

	static std::vector<Element*> deepCopyElements(Element* const*const mesh);


	unsigned _mesh_dimension;
	double _edge_length[2];
	std::string _name;
	std::vector<Node*> _nodes;
	std::vector<Element*> _elements;

}; /* class */

} /* namespace */

#endif /* MESH_H_ */

