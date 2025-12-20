// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "EdgeReturn.h"
#include "LineRule.h"

namespace MeshLib
{
/**
 * A 1d Edge or Line Element with 2 nodes.
 * \code
 *  0--------1
 * \endcode
 */
class LineRule2 : public LineRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 2u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::LINE2;

    /// Constant: Local node index table for edge
    constexpr static const unsigned edge_nodes[1][2] = {
        {0, 1}  // Edge 0
    };

    /// Edge rule
    using EdgeReturn = MeshLib::LinearEdgeReturn;
}; /* class */

}  // namespace MeshLib
