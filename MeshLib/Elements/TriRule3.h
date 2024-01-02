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
#include "TriRule.h"

namespace MeshLib
{

/**
 * This class represents a 2d triangle element with 3 nodes.
 *
 * The following sketch shows the node and edge numbering.
 * \anchor TriNodeAndEdgeNumbering
 * \code
 *
 *          2
 *          o
 *         / \
 *        /   \
 *      2/     \1
 *      /       \
 *     /         \
 *    0-----------1
 *          0
 *
 * \endcode
 */
class TriRule3 : public TriRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 3u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::TRI3;

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[3][2];

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::LinearEdgeReturn;

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[2])
    {
        return FaceRule::identifyFace<TriRule3>(element_nodes, nodes);
    }

}; /* class */

}  // namespace MeshLib
