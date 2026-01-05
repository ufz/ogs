// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "EdgeReturn.h"
#include "LineRule.h"

namespace MeshLib
{

/**
 * A 1d Edge or Line element with 3 nodes.
 * \code
 *  0----2----1
 * \endcode
 */
class LineRule3 : public LineRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 3u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::LINE3;

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[1][3] = {
        {0, 1, 2}  // Edge 0
    };

    /// Edge rule
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;
}; /* class */

}  // namespace MeshLib
