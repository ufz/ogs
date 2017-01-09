/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 2;

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[1][2];

    /// Edge rule
    typedef LinearEdgeReturn EdgeReturn;

    /**
     * \copydoc MeshLib::Element::isPntInElement()
     * @param nodes the nodes of the element.
     */
    static bool isPntInElement(Node const* const* nodes,
                               MathLib::Point3d const& pnt, double eps);

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
