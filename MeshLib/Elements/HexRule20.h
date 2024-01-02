/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    static const unsigned face_nodes[6][8];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[12][3];

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
