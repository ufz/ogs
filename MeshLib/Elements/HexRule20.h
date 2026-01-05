// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "EdgeReturn.h"
#include "HexRule.h"

namespace MeshLib
{

/**
 * A 20-nodes Hexahedron Element.
 * \code
 *
 *  Hex:
 *                14
 *          7-----------6
 *         /:          /|
 *        / :         / |
 *     15/  :        /13|
 *      / 19:       /   | 18
 *     /    : 12   /    |
 *    4-----------5     |
 *    |     :     | 10  |
 *    |     3.....|.....2
 *    |    .      |    /
 * 16 |   .       |17 /
 *    |11.        |  / 9
 *    | .         | /
 *    |.          |/
 *    0-----------1
 *          8
 *
 * \endcode
 */
class HexRule20 : public HexRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 20u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::HEX20;

    /// Constant: Local node index table for faces
    constexpr static const unsigned face_nodes[6][8] = {
        {0, 3, 2, 1, 11, 10, 9, 8},    // Face 0
        {0, 1, 5, 4, 8, 17, 12, 16},   // Face 1
        {1, 2, 6, 5, 9, 18, 13, 17},   // Face 2
        {2, 3, 7, 6, 10, 19, 14, 18},  // Face 3
        {3, 0, 4, 7, 11, 16, 15, 19},  // Face 4
        {4, 5, 6, 7, 12, 13, 14, 15}   // Face 5
    };

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[12][3] = {
        {0, 1, 8},   // Edge 0
        {1, 2, 9},   // Edge 1
        {2, 3, 10},  // Edge 2
        {0, 3, 11},  // Edge 3
        {4, 5, 12},  // Edge 4
        {5, 6, 13},  // Edge 5
        {6, 7, 14},  // Edge 6
        {4, 7, 15},  // Edge 7
        {0, 4, 16},  // Edge 8
        {1, 5, 17},  // Edge 9
        {2, 6, 18},  // Edge 10
        {3, 7, 19}   // Edge 11
    };

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return CellRule::identifyFace<HexRule20>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
