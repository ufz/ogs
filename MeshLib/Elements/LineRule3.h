/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshEnums.h"
#include "LineRule2.h"

namespace MeshLib
{

/**
 * A 1d Edge or Line element with 3 nodes.
 * @code
 *  0----2----1
 * @endcode
 */
class LineRule3 : public LineRule2
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 3u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::LINE3;

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[1][3];

    /// Edge rule
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;
}; /* class */

}  // namespace MeshLib
