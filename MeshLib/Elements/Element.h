/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Element class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <limits>
#include "MeshEnums.h"
#include "Mesh.h"
#include "MeshQuality/ElementErrorCode.h"

namespace MeshLib {

class Node;

/**
 * Virtual base class for mesh elements.
 */
class Element
{
	/* friend classes */
	friend class Mesh;//void Mesh::setElementInformationForNodes();

public:
	/// Compute the minimum and maximum squared edge length for this element
	virtual void computeSqrEdgeLengthRange(double &min, double &max) const;

	/**
	 * \brief Tries to add an element e as neighbour to this element.
	 * If the elements really are neighbours, the element is added to the
	 * neighbour-list and true is returned. Otherwise false is returned.
	 */
	bool addNeighbor(Element* e);

	/// Returns the length, area or volume of a 1D, 2D or 3D element
	virtual double getContent() const = 0;

	/**
	 * Get node with local index i where i should be at most the number
	 * of nodes of the element
	 * @param i local index of node, at most the number of nodes of the
	 * element that you can obtain with Element::getNNodes()
	 * @return a pointer to the appropriate (and constant, i.e. not
	 * modifiable by the user) instance of class Node or a NULL pointer
	 * @sa Element::getNodeIndex()
	 */
	const Node* getNode(unsigned i) const;

	/**
	 * (Re)Sets the node of the element.
	 * @param idx the index of the pointer to a node within the element
	 * @param node a pointer to a node
	 */
	void setNode(unsigned idx, Node* node);

	/// Get array of element nodes.
	Node* const* getNodes() const { return _nodes; }

	/// Get dimension of the mesh element.
	virtual unsigned getDimension() const = 0;

	/// Returns the i-th edge of the element.
	const Element* getEdge(unsigned i) const;

	/// Returns the i-th face of the element.
	virtual const Element* getFace(unsigned i) const = 0;

	/// Returns the ID of the element.
	virtual std::size_t getID() const { return this->_id; }

	/// Get the number of edges for this element.
	virtual unsigned getNEdges() const = 0;

	/// Get the number of nodes for face i.
	virtual unsigned getNFaceNodes(unsigned i) const = 0;

	/// Get the number of faces for this element.
	virtual unsigned getNFaces() const = 0;

	/// Get the specified neighbor.
	const Element* getNeighbor(unsigned i) const;

	/// Get the number of neighbors for this element.
	virtual unsigned getNNeighbors() const = 0;

	/**
	 * Returns the number of nodes. In dependency of the parameter
	 * the number of nodes for the geometric element is returned or
	 * the total number of nodes associated with this element
	 * is returned. The numbers can be different for instance if the
	 * element is used for higher order elements in finite element
	 * method.
	 * @param all (default = false)
	 * @return the number of nodes with respect to the parameter.
	 */
	virtual unsigned getNNodes(bool all = false) const = 0;

	/// Returns the position of the given node in the node array of this element.
	virtual unsigned getNodeIDinElement(const MeshLib::Node* node) const;

	/**
	 * Get the global index for the Node with local index i.
	 * The index i should be at most the number of nodes of the element.
	 * @param i local index of Node, at most the number of nodes of the
	 * element that you can obtain with Element::getNNodes()
	 * @return the global index or std::numeric_limits<unsigned>::max()
	 * @sa Element::getNode()
	 */
	unsigned getNodeIndex(unsigned i) const;

	/**
	 * Get the type of the mesh element in geometric context (as a MeshElemType-enum).
	 */
	virtual MeshElemType getGeomType() const = 0;

	/**
	 * Get the type of the element in context of the finite element method.
	 * @return a value of the enum FEMElemType::type
	 */
	virtual CellType getCellType() const = 0;

	/**
	 * Get the value for this element. The value can be used to store a link
	 * to external information (for instance an index of a field) like material groups.
	 */
	unsigned getValue() const { return _value; }

	/**
	 * Returns true if the element has zero length/area/volume.
	 */
	bool hasZeroVolume() const { return this->getContent() < std::numeric_limits<double>::epsilon(); }


	/// Returns true if these two indeces form an edge and false otherwise
	virtual bool isEdge(unsigned i, unsigned j) const = 0;

	/**
	 * Tests if the element is geometrically valid.
	 */
	virtual ElementErrorCode validate() const = 0;

	/**
	 * Set the index value for external information.
	 * @param value an unsigned value for linking with external information
	 */
	void setValue(unsigned value) { _value = value; }

	/// Returns true if elem is a neighbour of this element and false otherwise.
	bool hasNeighbor(Element* elem) const;

	/// Destructor
	virtual ~Element();

	/**
	 * Method clone is a pure virtual method in the abstract base class Element.
	 * It has to be implemented in the derived classes (for instance in class Hex).
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const = 0;

	/**
	 * Computes the length / area / volumen of this element. This is automatically
	 * done at initalisation time but can be repeated by calling this function at any time.
	 */
	virtual double computeVolume() = 0;


protected:
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Element(unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Return a specific edge node.
	virtual Node* getEdgeNode(unsigned edge_id, unsigned node_id) const = 0;

	/// Returns the ID of a face given an array of nodes.
	virtual unsigned identifyFace(Node* nodes[3]) const = 0;

	/// Sets the element ID.
	virtual void setID(std::size_t id) { this->_id = id; }


	Node** _nodes;
	std::size_t _id;
	/**
	 * this is an index for external additional information like materials
	 */
	unsigned _value;
	Element** _neighbors;
}; /* class */

} /* namespace */

#endif /* ELEMENT_H_ */

