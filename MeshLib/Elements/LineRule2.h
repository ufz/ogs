/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINERULE2_H_
#define LINERULE2_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "EdgeRule.h"
#include "EdgeReturn.h"

namespace MeshLib
{

/**
 * A 1d Edge or Line Element.
 * @code
 *  0--------1
 * @endcode
 */
class LineRule2 : public EdgeRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 2u;

    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 2u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::LINE;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::LINE2;

    /// Constant: The number of edges
    static const unsigned n_edges = 1;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 2;

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[1][2];

    /// Edge rule
    typedef NoEdgeReturn EdgeReturn;

    /**
     * Checks if a point is inside the element.
     * @param pnt a 3D MathLib::Point3d object
     * @param eps tolerance for numerical algorithm used or computing the property
     * @return true if the point is not outside the element, false otherwise
     */
    static bool isPntInElement(Node const* const* _nodes, MathLib::Point3d const& pnt, double eps);

    /**
     * Tests if the element is geometrically valid.
     */
    static ElementErrorCode validate(const Element* e);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const*, Node* nodes[1]);

    /// Calculates the length of a line
    static double computeVolume(Node const* const* _nodes);

}; /* class */

} /* namespace */

#endif /* LINERULE2_H_ */

