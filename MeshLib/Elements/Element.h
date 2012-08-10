/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Element.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include "MshEnums.h"
#include "Mesh.h"

namespace MeshLib {

class Node;

/**
 * Virtual base class for mesh elements.
 */
class Element
{
	/* friend functions: */
	friend class Mesh;//void Mesh::setElementInformationForNodes();
	//friend void Mesh::addElement(Element*);


public:
	/// Compute the minimum and maximum squared edge length for this element
	virtual void computeSqrEdgeLengthRange(double &min, double &max) const;

	/// Returns the length, area or volume of a 1D, 2D or 3D element
	virtual double getContent() const = 0;

	/// Get node with local index i.
	const Node* getNode(unsigned i) const;

	/**
	 * (Re)Sets the node of the element.
	 * @param idx the index of the pointer to a node within the element
	 * @param node a pointer to a node
	 */
	void setNode(unsigned idx, Node* node);

	/// Get array of element nodes.
	Node* const* getNodes() const { return _nodes; };

	/// Get dimension of the mesh element.
	virtual unsigned getDimension() const = 0;

	/// Returns the edge i of the element.
	const Element* getEdge(unsigned i) const;

	/// Returns the face i of the element.
	virtual const Element* getFace(unsigned i) const = 0;

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

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes() const = 0;

	/// Get the global index for the node with local index i.
	unsigned getNodeIndex(unsigned i) const;

	/// Get the type of the mesh element (as a MshElemType-enum).
	virtual MshElemType::type getType() const = 0;

	/// Get the value for this element.
	unsigned getValue() const { return _value; };

	bool hasNeighbor(Element* elem) const;

	/// Destructor
	virtual ~Element();

	/**
	 * Method clone is a pure virtual method in the abstract base class Element.
	 * It has to be implemented in the derived classes (for instance in class Hex).
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const = 0;

protected:
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Element(unsigned value = 0);

	/// Return a specific edge node.
	virtual Node* getEdgeNode(unsigned edge_id, unsigned node_id) const = 0;

	Node** _nodes;
	unsigned _value;
	Element** _neighbors;

private:

}; /* class */

} /* namespace */

#endif /* ELEMENT_H_ */

