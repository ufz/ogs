/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the TemplateTri class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATETRI_H_
#define TEMPLATETRI_H_

#include <array>

#include "GeoLib/AnalyticalGeometry.h"

#include "MathLib/MathTools.h"

#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"

#include "Face.h"
#include "Line.h"

namespace MeshLib {

namespace detail
{
class TriEdgeLinearNodes
{
protected:
	static constexpr unsigned _edge_nodes[3][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{0, 2}  // Edge 2
	};
};
} // end detail

/**
 * This class represents a 2d triangle element. The following sketch shows the node and edge numbering.
 * @anchor TriNodeAndEdgeNumbering
 * @code
 *
 *          2
 *          o
 *         / \
 *        /   \
 *      2/     \1
 *      /       \
 *     /         \
 *    0-----------1
 *          0
 *
 * @endcode
 */
template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES = detail::TriEdgeLinearNodes>
class TemplateTri : public Face, public EDGENODES
{
protected:
	using EDGENODES::_edge_nodes;

public:
	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = NNODES;

	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes = 3u;

	/// Constructor with an array of mesh nodes.
	TemplateTri(Node* nodes[NNODES], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs a triangle from array of Node pointers.
	TemplateTri(std::array<Node*, NNODES> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Copy constructor
	TemplateTri(const TemplateTri &tri);

	/// Destructor
	virtual ~TemplateTri();

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 3; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 3; };

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
	 * Method returns the type of the element. In this case TRIANGLE will be returned.
	 * @return MeshElemType::TRIANGLE
	 */
	virtual MeshElemType getGeomType() const { return MeshElemType::TRIANGLE; }

	/**
	 * Get the type of the element in context of the finite element method.
	 * @return a value of the enum CellType
	 */
	virtual CellType getCellType() const { return CELLTRITYPE; }

	/// Returns true if these two indices form an edge and false otherwise
	bool isEdge(unsigned idx1, unsigned idx2) const;

	/**
	 * Checks if a point is inside the element.
	 * @param pnt a 3D GeoLib::Point object
	 * @param eps tolerance for numerical algorithm used or computing the property
	 * @return true if the point is not outside the element, false otherwise
	 */
	bool isPntInElement(GeoLib::Point const& pnt, double eps = std::numeric_limits<double>::epsilon()) const;

	/**
	 * Tests if the element is geometrically valid
	 * @param check_zero_volume indicates if area == 0 should be checked
	 * @return error code (0 = okay, 1 = zero volume)
	 */
	virtual ElementErrorCode validate() const;


	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateTri instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const
	{
		return new TemplateTri<NNODES,CELLTRITYPE>(*this);
	}

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const;

protected:
	/// Calculates the area of the triangle by returning half of the area of the corresponding parallelogram.
	double computeVolume()
	{
		return GeoLib::calcTriangleArea(*_nodes[0], *_nodes[1], *_nodes[2]);
	}

protected:
	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const
	{
		return _nodes[_edge_nodes[edge_id][node_id]];
	}

}; /* class */

} /* namespace */

#include "TemplateTri-impl.h"

#endif /* TEMPLATETRI_H_ */

