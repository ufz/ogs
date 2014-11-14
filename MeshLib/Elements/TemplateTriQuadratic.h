/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATETRIQUADRATIC_H_
#define TEMPLATETRIQUADRATIC_H_

#include <array>
#include "MeshLib/MeshEnums.h"
#include "TemplateTri.h"

namespace MeshLib {

namespace detail
{
class TriEdgeQuadraticNodes
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
template <unsigned NNODES, CellType CELLQUADTYPE>
class TemplateTriQuadratic : public TemplateTri<NNODES, CELLQUADTYPE, detail::TriEdgeQuadraticNodes>
{
public:
	/// Constructor with an array of mesh nodes.
	TemplateTriQuadratic(Node* nodes[NNODES], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs an edge from array of Node pointers.
	TemplateTriQuadratic(std::array<Node*, NNODES> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs a quad from NNODES of Nodes initializing Face with
	/// value = 0.
	TemplateTriQuadratic(Node* n0, Node* n1, Node* n2, Node* n3, ...);

	/// Copy constructor
	TemplateTriQuadratic(const TemplateTriQuadratic &tri);

	/// Destructor
	virtual ~TemplateTriQuadratic() {}

	/// Returns the i-th edge of the element.
	virtual const Element* getEdge(unsigned i) const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateQuad instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const;

}; /* class */

} /* namespace */

#include "TemplateTriQuadratic-impl.h"

#endif /* TEMPLATETRIQUADRATIC_H_ */

