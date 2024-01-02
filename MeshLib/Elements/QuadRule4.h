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
#include "QuadRule.h"

namespace MeshLib
{

/**
 * This class represents a 2d quadrilateral element with 4 nodes.
 * The following sketch shows the node and edge numbering.
 * \anchor QuadNodeAndEdgeNumbering
 * \code
 *              2
 *        3-----------2
 *        |           |
 *        |           |
 *       3|           |1
 *        |           |
 *        |           |
 *        0-----------1
 *              0
 * \endcode
 */
class QuadRule4 : public QuadRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 4u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::QUAD4;

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[4][2];

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::LinearEdgeReturn;

    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return FaceRule::identifyFace<QuadRule4>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
