/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEELEMENT_H_
#define TEMPLATEELEMENT_H_

#include <array>
#include <limits>

#include "MathLib/Point3d.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshQuality/ElementErrorCode.h"

namespace MeshLib
{

/**
 * Template for implementing mesh element classes
 *
 * \tparam T_BASE         Base element class, e.g. Face, Cell
 * \tparam ELEMENT_RULE   Geometrical and topological rules of the element
 */
template <class T_BASE, class ELEMENT_RULE>
class TemplateElement : public T_BASE
{
public:
	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = ELEMENT_RULE::n_all_nodes;

	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes = ELEMENT_RULE::n_base_nodes;

	/// Constant: The dimension of this element
	static const unsigned dimension = ELEMENT_RULE::dimension;

	/**
	 * Constructor with an array of mesh nodes.
	 *
	 * @param nodes  an array of pointers of mesh nodes which form this element
	 * @param value  element value, e.g. material ID
	 * @param id     element id
	 */
	TemplateElement(Node* nodes[n_all_nodes], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/**
	 * Constructor with an array of mesh nodes
	 *
	 * @param nodes  an array of pointers of mesh nodes which form this element
	 * @param value  element value, e.g. material ID
	 * @param id     element id
	 */
	TemplateElement(std::array<Node*, n_all_nodes> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Copy constructor
	TemplateElement(const TemplateElement &e);

	/// Destructor
	virtual ~TemplateElement() {}

	/// Returns a copy of this object.
	virtual Element* clone() const
	{
		return new TemplateElement(*this);
	}

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return dimension; }

	/// Returns the edge i of the element.
	const Element* getEdge(unsigned i) const { return ELEMENT_RULE::EdgeReturn::getEdge(this, i); }

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const { return ELEMENT_RULE::getFace(this, i); }

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return ELEMENT_RULE::n_edges; }

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { return ELEMENT_RULE::getNFaceNodes(i); }

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return ELEMENT_RULE::n_faces; }

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return ELEMENT_RULE::n_neighbors; }

	/// Get the number of linear nodes for this element.
	virtual unsigned getNBaseNodes() const { return n_base_nodes; }

	/// Get the number of all nodes for this element.
	virtual unsigned getNNodes() const { return n_all_nodes; }

	/// Get the type of this element.
	virtual MeshElemType getGeomType() const { return ELEMENT_RULE::mesh_elem_type; }

	/// Get the FEM type of this element.
	virtual CellType getCellType() const { return ELEMENT_RULE::cell_type; }

	/// Returns true if these two indices form an edge and false otherwise
	bool isEdge(unsigned idx1, unsigned idx2) const;

	/**
	 * Checks if a point is inside the element.
	 * @param pnt a 3D GeoLib::Point object
	 * @param eps tolerance for numerical algorithm used or computing the property
	 * @return true if the point is not outside the element, false otherwise
	 */
	bool isPntInElement(MathLib::Point3d const& pnt, double eps = std::numeric_limits<double>::epsilon()) const
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
		if (getNEdges()>0)
			return const_cast<Node*>(this->_nodes[ELEMENT_RULE::edge_nodes[edge_id][node_id]]);
		else
			return nullptr;
	}

};

} // MeshLib

#include "TemplateElement-impl.h"

#endif /* TEMPLATEELEMENT_H_ */

