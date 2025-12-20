// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "EdgeReturn.h"
#include "Element.h"
#include "TetRule.h"

namespace MeshLib
{

/**
 * This class represents a 3d tetrahedron element with 4 nodes.
 * The following sketch shows the node and edge numbering.
 * \anchor Tetrahedron4NodeAndEdgeNumbering
 * \code
 *          3
 *         /|\
 *        / | \
 *      3/  |  \5
 *      /   |4  \
 *     /    |    \
 *    0.....|.....2
 *     \    |  2 /
 *      \   |   /
 *      0\  |  /1
 *        \ | /
 *         \|/
 *          1
 *
 * \endcode
 */
class TetRule4 : public TetRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 4u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::TET4;

    /// Constant: Local node index table for faces
    constexpr static const unsigned face_nodes[4][3] = {
        {0, 2, 1},  // Face 0
        {0, 1, 3},  // Face 1
        {1, 2, 3},  // Face 2
        {2, 0, 3}   // Face 3
    };

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[6][2] = {
        {0, 1},  // Edge 0
        {1, 2},  // Edge 1
        {0, 2},  // Edge 2
        {0, 3},  // Edge 3
        {1, 3},  // Edge 4
        {2, 3}   // Edge 5
    };

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::LinearEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return CellRule::identifyFace<TetRule4>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
