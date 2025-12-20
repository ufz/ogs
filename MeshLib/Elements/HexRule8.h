// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "EdgeReturn.h"
#include "HexRule.h"

namespace MeshLib
{
/**
 * A 8-nodes Hexahedron Element.
 * \code
 *
 *  Hex:
 *                6
 *          7-----------6
 *         /:          /|
 *        / :         / |
 *      7/  :        /5 |
 *      / 11:       /   | 10
 *     /    : 4    /    |
 *    4-----------5     |
 *    |     :     | 2   |
 *    |     3.....|.....2
 *    |    .      |    /
 *  8 |   .       |9  /
 *    | 3.        |  / 1
 *    | .         | /
 *    |.          |/
 *    0-----------1
 *          0
 *
 * \endcode
 */
class HexRule8 : public HexRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 8u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::HEX8;

    /// Constant: Local node index table for faces
    constexpr static const unsigned face_nodes[6][4] = {
        {0, 3, 2, 1},  // Face 0
        {0, 1, 5, 4},  // Face 1
        {1, 2, 6, 5},  // Face 2
        {2, 3, 7, 6},  // Face 3
        {3, 0, 4, 7},  // Face 4
        {4, 5, 6, 7}   // Face 5
    };

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[12][2] = {
        {0, 1},  // Edge 0
        {1, 2},  // Edge 1
        {2, 3},  // Edge 2
        {0, 3},  // Edge 3
        {4, 5},  // Edge 4
        {5, 6},  // Edge 5
        {6, 7},  // Edge 6
        {4, 7},  // Edge 7
        {0, 4},  // Edge 8
        {1, 5},  // Edge 9
        {2, 6},  // Edge 10
        {3, 7}   // Edge 11
    };

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::LinearEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return CellRule::identifyFace<HexRule8>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
