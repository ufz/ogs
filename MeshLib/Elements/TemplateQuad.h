/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the TemplateQuad class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEQUAD_H_
#define TEMPLATEQUAD_H_

#include <array>
#include "MeshEnums.h"
#include "Face.h"

namespace MeshLib {

/**
 * This class represents a 2d quadliteral element. The following sketch shows the node and edge numbering.
 * @anchor QuadNodeAndEdgeNumbering
 * @code
 *              2
 *        3-----------2
 *        |           |
 *        |           |
 *       3|           |1
 *        |           |
 *        |           |
 *        0-----------1
 *              0
 * @endcode
 */
template <unsigned NNODES, CellType CELLQUADTYPE>
class TemplateQuad : public Face
{
public:
	/// Constant: The number of all nodes for this element
	static const unsigned _n_all_nodes;

	/// Constant: The number of base nodes for this element
	static const unsigned _n_base_nodes;

	/// Constructor with an array of mesh nodes.
	TemplateQuad(Node* nodes[NNODES], unsigned value = 0);

	/// Constructs an edge from array of Node pointers.
	TemplateQuad(std::array<Node*, NNODES> const& nodes, unsigned value = 0);

	/// Constructs a quad from NNODES of Nodes initializing Face with
	//  value = 0.
	TemplateQuad(Node* n0, Node* n1, Node* n2, Node* n3, ...);

	/// Copy constructor
	TemplateQuad(const TemplateQuad &quad);

	/// Destructor
	virtual ~TemplateQuad();

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 4; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 4; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(bool all = false) const
	{
		return all ? _n_all_nodes : _n_base_nodes;
	}

	/**
	 * Method returns the type of the element. In this case QUAD will be returned.
	 * @return MeshElemType::QUAD
	 */
	virtual MeshElemType getGeomType() const { return MeshElemType::QUAD; }

	/**
	 * Get the type of the element in context of the finite element method.
	 * @return a value of the enum CellType
	 */
	virtual CellType getCellType() const { return CELLQUADTYPE; }

	/// Returns true if these two indeces form an edge and false otherwise
	bool isEdge(unsigned i, unsigned j) const;

	/**
	 * Check if the 3d GeoLib::Point is inside of the quad element.
	 * @param pnt the 3d GeoLib::Point object
	 * @param eps tolerance for numerical algorithm used or computing the property
	 * @return true if the point is inside the element, false otherwise
	 */
	virtual bool isPntInside(GeoLib::Point const& pnt, double eps = std::numeric_limits<double>::epsilon()) const;


	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateQuad instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

	/**
	 * This method should be called after at least two nodes of the quad
	 * element are collapsed. As a consequence of the node collapsing an edge
	 * of the quad will be collapsed. If one of the edges (see
	 * sketch @ref PyramidNodeAndEdgeNumbering) is collapsed we obtain a
	 * triangle. In this case the method will create the appropriate
	 * object of class Tri.
	 * @return a Tri object or NULL
	 */
	virtual Element* reviseElement() const;

protected:
	/// Calculates the area of a convex quadliteral by dividing it into two triangles.
	double computeVolume();

protected:
	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

	static const unsigned _edge_nodes[4][2];
}; /* class */

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad<NNODES, CELLQUADTYPE>::_edge_nodes[4][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}  // Edge 3
};

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad<NNODES, CELLQUADTYPE>::_n_all_nodes = NNODES;

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad<NNODES, CELLQUADTYPE>::_n_base_nodes = 4;

} /* namespace */

#include "TemplateQuad.tpp"

#endif /* TEMPLATEQUAD_H_ */

