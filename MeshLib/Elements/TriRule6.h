/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshEnums.h"
#include "TriRule3.h"
#include "EdgeReturn.h"

namespace MeshLib
{

/**
 * This class represents a 2d triangle element with 6 nodes.
 *
 * The following sketch shows the node and edge numbering.
 * @anchor Tri6NodeAndEdgeNumbering
 * @code
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
 * @endcode
 */
class TriRule6 : public TriRule3
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

}; /* class */

} /* namespace */
