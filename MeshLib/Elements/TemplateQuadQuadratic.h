/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEQUADQUADRATIC_H_
#define TEMPLATEQUADQUADRATIC_H_

#include <array>
#include "MeshLib/MeshEnums.h"
#include "TemplateQuad.h"

namespace MeshLib {

namespace detail
{
class QuadEdgeQuadraticNodes
{
protected:
	static constexpr unsigned _edge_nodes[4][3] = {
		{0, 1, 4}, // Edge 0
		{1, 2, 5}, // Edge 1
		{2, 3, 6}, // Edge 2
		{0, 3, 7}  // Edge 3
	};
};
} // end detail

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
class TemplateQuadQuadratic : public TemplateQuad<NNODES, CELLQUADTYPE, detail::QuadEdgeQuadraticNodes>
{
public:
	/// Constructor with an array of mesh nodes.
	TemplateQuadQuadratic(Node* nodes[NNODES], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs an edge from array of Node pointers.
	TemplateQuadQuadratic(std::array<Node*, NNODES> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs a quad from NNODES of Nodes initializing Face with
	//  value = 0.
	TemplateQuadQuadratic(Node* n0, Node* n1, Node* n2, Node* n3, ...);

	/// Copy constructor
	TemplateQuadQuadratic(const TemplateQuadQuadratic &quad);

	/// Destructor
	virtual ~TemplateQuadQuadratic() {}

	/// Returns the i-th edge of the element.
	virtual const Element* getEdge(unsigned i) const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateQuad instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

}; /* class */

} /* namespace */

#include "TemplateQuadQuadratic-impl.h"

#endif /* TEMPLATEQUADQUADRATIC_H_ */

