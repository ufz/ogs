/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEELEMENT_H_
#define TEMPLATEELEMENT_H_

#include "Element.h"

namespace MeshLib {

template <class T_BASE, class ELEMENT_RULE>
class TemplateElement : public T_BASE
{
public:
//	/// Constant: Dimension of this mesh element
//	static const unsigned dimension = ELEMENT_RULE::dimension;

	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = ELEMENT_RULE::n_all_nodes;

	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes = ELEMENT_RULE::n_base_nodes;

//	/// Returns the length, area or volume of a 1D, 2D or 3D element
//	double getContent() const { return _content; };

//	/// Get dimension of the mesh element.
//	unsigned getDimension() const { return this->dimension; }

//	/// Get the volume of this 3d element.
//	virtual double getVolume() const { return _content; };

	/// Constructor with an array of mesh nodes.
	TemplateElement(Node* nodes[n_all_nodes], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max())
	: T_BASE(value, id)
	{
		this->_nodes = nodes;
		this->_neighbors = new Element*[getNNeighbors()];
		std::fill(this->_neighbors, this->_neighbors + getNNeighbors(), nullptr);
		this->_volume = ELEMENT_RULE::computeVolume(this->_nodes);
	}

	/// Constructs a hex from array of Node pointers.
	TemplateElement(std::array<Node*, n_all_nodes> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max())
	: T_BASE(value, id)
	{
		this->_nodes = new Node*[n_all_nodes];
		std::copy(nodes.begin(), nodes.end(), this->_nodes);
		this->_neighbors = new Element*[getNNeighbors()];
		std::fill(this->_neighbors, this->_neighbors + getNNeighbors(), nullptr);
		this->_volume = ELEMENT_RULE::computeVolume(this->_nodes);
	}

	/// Copy constructor
	TemplateElement(const TemplateElement &e)
	: T_BASE(e.getValue(), e.getID())
	{
		this->_nodes = new Node*[n_all_nodes];
		for (unsigned i=0; i<n_all_nodes; i++)
			this->_nodes[i] = e._nodes[i];
		this->_neighbors = new Element*[getNNeighbors()];
		for (unsigned i=0; i<getNNeighbors(); i++)
			this->_neighbors[i] = e._neighbors[i];
		this->_volume = e.getContent();
	}

	/// Destructor
	virtual ~TemplateElement() {}

	/**
	 * This method is pure virtual and is inherited from class @sa Element.
	 * It has to be implemented in the derived classes of class Cell!
	 * @return a copy of the object
	 */
	virtual Element* clone() const
	{
		return new TemplateElement(*this);
	}

//	/**
//	 * Checks if the node order of an element is correct by testing surface normals.
//	 * For 3D elements true is returned if the normals of all faces points away from the centre of
//	 * the element.
//	 * Note: This method might give wrong results if something else is wrong with the element
//	 * (non-planar faces, non-convex geometry, possibly zero volume) which causes the calculated
//	 * center of gravity to lie outside of the actual element
//	 */
//	virtual bool testElementNodeOrder() const { return ELEMENT_RULE::testElementNodeOrder(this); }

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const { return ELEMENT_RULE::getFace(this->_nodes, i); }

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return ELEMENT_RULE::n_edges; }

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { return ELEMENT_RULE::getNFaceNodes(i); }

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return ELEMENT_RULE::n_faces; }

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return ELEMENT_RULE::n_neighbors; }

	/// Get the number of linear nodes for this element.
	virtual unsigned getNBaseNodes() const
	{
		return n_base_nodes;
	}

	/// Get the number of all nodes for this element.
	virtual unsigned getNNodes() const
	{
		return n_all_nodes;
	}

	/**
	 * Method returns the type of the element. In this case HEXAHEDRON will be returned.
	 * @return MeshElemType::HEXAHEDRON
	 */
	virtual MeshElemType getGeomType() const { return ELEMENT_RULE::mesh_elem_type; }

	/**
	 * Method returns the FEM type of the element.
	 * @return
	 */
	virtual CellType getCellType() const { return ELEMENT_RULE::cell_type; }

	/// Returns true if these two indices form an edge and false otherwise
	bool isEdge(unsigned idx1, unsigned idx2) const
	{
		for (unsigned i(0); i<getNEdges(); i++)
		{
			if (ELEMENT_RULE::_edge_nodes[i][0]==idx1 && ELEMENT_RULE::_edge_nodes[i][1]==idx2) return true;
			if (ELEMENT_RULE::_edge_nodes[i][1]==idx1 && ELEMENT_RULE::_edge_nodes[i][0]==idx2) return true;
		}
		return false;
	}

	/**
	 * Checks if a point is inside the element.
	 * @param pnt a 3D GeoLib::Point object
	 * @param eps tolerance for numerical algorithm used or computing the property
	 * @return true if the point is not outside the element, false otherwise
	 */
	bool isPntInElement(GeoLib::Point const& pnt, double eps = std::numeric_limits<double>::epsilon()) const
	{
		return ELEMENT_RULE::isPntInElement(this->_nodes, pnt, eps);
	}

	/**
	 * Tests if the element is geometrically valid.
	 * @param check_zero_volume indicates if volume == 0 should be checked
	 */
	virtual ElementErrorCode validate() const
	{
		return ELEMENT_RULE::validate(this);
	}

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const
	{
		return ELEMENT_RULE::identifyFace(this->_nodes, nodes);
	};

	/// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
	virtual double computeVolume() {return ELEMENT_RULE::computeVolume(this->_nodes);}

	/// Return a specific edge node.
	virtual inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const
	{
		return const_cast<Node*>(this->_nodes[ELEMENT_RULE::_edge_nodes[edge_id][node_id]]);
	}


protected:
//	double _content;

}; /* class */

}

#endif /* TEMPLATEELEMENT_H_ */

