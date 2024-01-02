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
#include "Element.h"
#include "PrismRule.h"

namespace MeshLib
{

/**
 * This class represents a 3d prism element with 6 nodes.
 * The following sketch shows the node and edge numbering.
 * \anchor Prism6NodeAndEdgeNumbering
 * \code
 *            5
 *           / \
 *          / : \
 *        8/  :  \7
 *        /   :5  \
 *       /    :  6 \
 *      3-----------4
 *      |     :     |
 *      |     2     |
 *      |    . .    |
 *     3|   .   .   |4
 *      | 2.     .1 |
 *      | .       . |
 *      |.         .|
 *      0-----------1
 *            0
 *
 * \endcode
 */
class PrismRule6 : public PrismRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 6u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::PRISM6;

    /// Constant: Local node index table for faces
    static const unsigned face_nodes[5][4];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[9][2];

    /// Constant: Table for the number of nodes for each face
    static const unsigned n_face_nodes[5];

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::LinearEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return CellRule::identifyFace<PrismRule6>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
