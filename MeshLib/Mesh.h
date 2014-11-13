/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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
#include <map>

#include <boost/any.hpp>
#include <boost/optional.hpp>

#include "logog/include/logog.hpp"

#include "BaseLib/Counter.h"

#include "MeshEnums.h"
#include "Location.h"

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
	/// @param name          Mesh name.
	/// @param nodes         A vector of mesh nodes. In case nonlinear nodes are involved, one should
	///                      put them after line ones in the vector and set "n_base_nodes" argument.
	/// @param elements      An array of mesh elements.
	/// @param n_base_nodes  The number of base nodes. This is an optional parameter for nonlinear case.
	///                      If the parameter is set to zero, we consider there are no nonlinear nodes.
	Mesh(const std::string &name,
	     const std::vector<Node*> &nodes,
	     const std::vector<Element*> &elements,
	     const std::size_t n_base_nodes = 0);

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
	double getMinEdgeLength() const { return _edge_length.first; }

	/// Get the maximum edge length over all elements of the mesh.
	double getMaxEdgeLength() const { return _edge_length.second; }

	/// Get the minimum node distance in the mesh.
    /// The value is calculated from element-wise minimum node distances.
	double getMinNodeDistance() const { return _node_distance.first; }

	/// Get the maximum node distance over all elements of the mesh.
    /// The value is calculated from element-wise maximum node distances.
	double getMaxNodeDistance() const { return _node_distance.second; }

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

	/// Get the number of base nodes
	std::size_t getNBaseNodes() const { return _n_base_nodes; }

	/// Return true if the given node is a basic one (i.e. linear order node)
	bool isBaseNode(std::size_t node_idx) const {return node_idx < _n_base_nodes; }

	/// Return true if the mesh has any nonlinear nodes
	bool isNonlinear() const { return (getNNodes() != getNBaseNodes()); }

	/// Method to get a vector of property values.
	template <typename T>
	boost::optional<std::vector<T> const&>
	getProperty(std::string const& name) const
	{
		std::map<std::string, boost::any>::const_iterator it(
			_properties.find(name)
		);
		if (it != _properties.end()) {
			try {
				boost::any_cast<std::vector<T> const&>(it->second);
				return boost::any_cast<std::vector<T> const&>(it->second);
			} catch (boost::bad_any_cast const&) {
				ERR("A property with the desired data type is not available.");
				return boost::optional<std::vector<T> const&>();
			}
		} else {
			return boost::optional<std::vector<T> const&>();
		}
	}

	/// Method to store a vector of property values assigned to a property name.
	/// Since the implementation makes no assumption about the number of data
	/// items stored within the vector, it is possible either to use a small
	/// number of properties where each particular property can be assigned to
	/// several mesh items. In contrast to this it is possible to have a
	/// separate value for each mesh item.
	/// The user has to ensure the correct usage of the vector later on.
	template <typename T>
	void addProperty(std::string const& name, std::vector<T> const& property)
	{
		std::map<std::string, boost::any>::const_iterator it(
			_properties.find(name)
		);
		if (it != _properties.end()) {
			WARN("A property of the name \"%s\" already assigned to the mesh.",
				name.c_str());
			return;
		}
		_properties[name] = boost::any(property);
	}

protected:
	/// Set the minimum and maximum length over the edges of the mesh.
	void calcEdgeLengthRange();
	/// Set the minimum and maximum node distances within elements.
	void calcNodeDistanceRange();

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
	/// Note: Using this implementation, an element e can only have neighbors that have the same dimensionality as e.
	void setElementNeighbors();

	void setNodesConnectedByEdges();

	void setNodesConnectedByElements();

	std::size_t const _id;
	unsigned _mesh_dimension;
	/// The minimal and maximal edge length over all elements in the mesh
	std::pair<double, double> _edge_length;
	/// The minimal and maximal distance of nodes within an element over all elements in the mesh
	std::pair<double, double> _node_distance;
	std::string _name;
	std::vector<Node*> _nodes;
	std::vector<Element*> _elements;
	std::size_t _n_base_nodes;

	struct PropertyKeyType
	{
		PropertyKeyType(std::string const& n, MeshItemType t)
			: name(n), mesh_item_type(t)
		{}

		std::string name;
		MeshItemType mesh_item_type;

		bool operator<(PropertyKeyType const& other) const
		{
			if (name.compare(other.name) == 0) {
				return mesh_item_type < other.mesh_item_type;
			}
			return name.compare(other.name) < 0 ? true : false;
		}
	};

	/// A mapping from property's name to the stored object of any type.
	/// See addProperty() and getProperty() documentation.
	std::map<PropertyKeyType, boost::any> _properties;

}; /* class */

} /* namespace */

#endif /* MESH_H_ */

