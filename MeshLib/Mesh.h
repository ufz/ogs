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

#include <boost/optional.hpp>

#include "MshEnums.h"

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
	 * This should have been previously calculated using the Element::computeSqrEdgeLengthRange(min, max)
	 * function or by some other means.
	 */
	void setEdgeLengthRange(const double &min_length, const double &max_length);

	/// Changes the name of the mesh.
	void setName(const std::string &name) { this->_name = name; };

	/// Get id of the mesh
	std::size_t getID() const {return _id; }

	/**
	 * Adds for each element a property associated with prop_name. The
	 * properties are stored in the vector elem_props. The number of entries
	 * of this vector should be the same as the number of elements for scalar
	 * values and a multiple of the number of elements for vector or matrix
	 * properties. This method adds properties that can be expressed deploying
	 * floating point numbers.
	 * @param prop_name The name of the property, for instance permeability.
	 * @param elem_props The vector containing the property values.
	 */
	void addPropertyVec(std::string const& prop_name, std::vector<double> const& elem_props);

	/**
	 * Adds for each element a property associated with prop_name. The
	 * properties are stored in the vector elem_props. The number of entries
	 * of this vector should be the same as the number of elements for scalar
	 * values and a multiple of the number of elements for vector or matrix
	 * properties. This method adds properties that can be expressed deploying
	 * integers.
	 * @param prop_name The name of the property, for instance MateriaID.
	 * @param elem_props The vector containing the property values.
	 */
	void addPropertyVec(std::string const& prop_name, std::vector<unsigned> const& elem_props);

	/**
	 * Get the vector of properties associated with the name prop_name.
	 * If there is not a vector associate with the given name an
	 * invalid_argument exception is thrown.
	 * @param prop_name The name of the property.
	 * @return Vector containing the properties.
	 */
	boost::optional<std::vector<double> const&>
	getDoublePropertyVec(std::string const& prop_name) const;

	/**
	 * Get the vector of properties associated with the name prop_name.
	 * If there is not a vector associate with the given name an
	 * invalid_argument exception is thrown.
	 * @param prop_name The name of the property.
	 * @return Vector containing the properties.
	 */
	boost::optional<std::vector<unsigned> const&>
	getUnsignedPropertyVec(std::string const& prop_name) const;

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

	std::size_t const _id;
	unsigned _mesh_dimension;
	double _edge_length[2];
	std::string _name;
	std::vector<Node*> _nodes;
	std::vector<Element*> _elements;

	std::vector<std::pair<std::string, std::vector<double> > > _double_prop_vecs;
	std::vector<std::pair<std::string, std::vector<unsigned> > > _unsigned_prop_vecs;
}; /* class */

} /* namespace */

#endif /* MESH_H_ */

