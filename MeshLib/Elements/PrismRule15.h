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
#include "PrismRule.h"

namespace MeshLib
{

/**
 * This class represents a 3d prism element with 15 nodes. The following sketch
 * shows the node and edge numbering.
 * \anchor PrismNodeAndEdgeNumbering
 * \code
 *            5
 *           / \
 *          / : \
 *       11/  :  \10
 *        /   :14 \
 *       /    :    \
 *      3------9----4
 *      |     :     |
 *      |     2     |
 *      |    . .    |
 *    12|   .   .   |13
 *      | 8.     .7 |
 *      | .       . |
 *      |.         .|
 *      0-----------1
 *            6
 *
 * \endcode
 */
class PrismRule15 : public PrismRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 15u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::PRISM15;

    /// Constant: Local node index table for faces
    constexpr static const unsigned face_nodes[5][8] = {
        {0, 2, 1, 8, 7, 6, 99, 99},   // Face 0
        {0, 1, 4, 3, 6, 13, 9, 12},   // Face 1
        {1, 2, 5, 4, 7, 14, 10, 13},  // Face 2
        {2, 0, 3, 5, 8, 12, 11, 14},  // Face 3
        {3, 4, 5, 9, 10, 11, 99, 99}  // Face 4
    };

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[9][3] = {
        {0, 1, 6},   // Edge 0
        {1, 2, 7},   // Edge 1
        {0, 2, 8},   // Edge 2
        {0, 3, 12},  // Edge 3
        {1, 4, 13},  // Edge 4
        {2, 5, 14},  // Edge 5
        {3, 4, 9},   // Edge 6
        {4, 5, 10},  // Edge 7
        {3, 5, 11}   // Edge 8
    };

    /// Constant: Table for the number of nodes for each face
    static const unsigned n_face_nodes[5];

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return CellRule::identifyFace<PrismRule15>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
