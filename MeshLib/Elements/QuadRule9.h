// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "EdgeReturn.h"
#include "MeshLib/MeshEnums.h"
#include "QuadRule.h"

namespace MeshLib
{

/**
 * This class represents a 2d quadrilateral element with 9 nodes.
 * The following sketch shows the node and edge numbering.
 * \anchor Quad9NodeAndEdgeNumbering
 * \code
 *              2
 *        3-----6-----2
 *        |           |
 *        |           |
 *      3 7     8     5 1
 *        |           |
 *        |           |
 *        0-----4-----1
 *              0
 * \endcode
 */
class QuadRule9 : public QuadRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 9u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::QUAD9;

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[4][3] = {
        {0, 1, 4},  // Edge 0
        {1, 2, 5},  // Edge 1
        {2, 3, 6},  // Edge 2
        {3, 0, 7}   // Edge 3
    };

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return FaceRule::identifyFace<QuadRule9>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
