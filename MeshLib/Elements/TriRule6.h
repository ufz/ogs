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
#include "MeshLib/MeshEnums.h"
#include "TriRule.h"

namespace MeshLib
{

/**
 * This class represents a 2d triangle element with 6 nodes.
 *
 * The following sketch shows the node and edge numbering.
 * \anchor Tri6NodeAndEdgeNumbering
 * \code
 *
 *          2
 *          o
 *         / \
 *        /   \
 *      5/     \4
 *      /       \
 *     /         \
 *    0-----------1
 *          3
 *
 * \endcode
 */
class TriRule6 : public TriRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 6u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::TRI6;

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[3][3];

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[2])
    {
        return FaceRule::identifyFace<TriRule6>(element_nodes, nodes);
    }

}; /* class */

}  // namespace MeshLib
