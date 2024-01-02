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
    static const unsigned edge_nodes[1][2];

    /// Edge rule
    using EdgeReturn = MeshLib::LinearEdgeReturn;
}; /* class */

}  // namespace MeshLib
