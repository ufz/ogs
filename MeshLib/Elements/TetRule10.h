/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "EdgeReturn.h"
#include "Element.h"
#include "TetRule.h"

namespace MeshLib
{

/**
 * This class represents a 3d tetrahedron element with 10 nodes. The following
 * sketch shows the node and edge numbering.
 * \anchor TetrahedronNodeAndEdgeNumbering
 * \code
 *          3
 *         /|\
 *        / | \
 *      7/  |  \9
 *      /   |8  \
 *     /    |    \
 *    0.....|.....2
 *     \    |  2 /
 *      \   |   /
 *      4\  |  /5
 *        \ | /
 *         \|/
 *          1
 *
 * \endcode
 */
class TetRule10 : public TetRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 10u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::TET10;

    /// Constant: Local node index table for faces
    constexpr static const unsigned face_nodes[4][6] = {
        {0, 2, 1, 6, 5, 4},  // Face 0
        {0, 1, 3, 4, 8, 7},  // Face 1
        {1, 2, 3, 5, 9, 8},  // Face 2
        {2, 0, 3, 6, 7, 9}   // Face 3
    };

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[6][3] = {
        {0, 1, 4},  // Edge 0
        {1, 2, 5},  // Edge 1
        {0, 2, 6},  // Edge 2
        {0, 3, 7},  // Edge 3
        {1, 3, 8},  // Edge 4
        {2, 3, 9}   // Edge 5
    };

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return CellRule::identifyFace<TetRule10>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
