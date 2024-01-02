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
    static const unsigned face_nodes[6][4];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[12][2];

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
