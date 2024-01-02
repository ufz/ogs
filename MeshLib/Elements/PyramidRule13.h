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
#include "PyramidRule.h"

namespace MeshLib
{

/**
 * This class represents a 3d pyramid element with 13 nodes.
 * The following sketch shows the node numbering.
 * \anchor Pyramid13NodeNumbering
 * \code
 *
 *               4
 *             //|\
 *            // | \
 *        12 //  |  \11
 *          //   |10 \
 *         //    |    \
 *        3/.... |.....2
 *       ./      |  7 /
 *      ./9      |   /
 *    8./        |  /6
 *    ./         | /
 *   ./          |/
 *  0------------1
 *        5
 * \endcode

 */
class PyramidRule13 : public PyramidRule
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 13u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::PYRAMID13;

    /// Constant: Local node index table for faces
    static const unsigned face_nodes[5][8];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[8][3];

    /// Constant: Table for the number of nodes for each face
    static const unsigned n_face_nodes[5];

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[3])
    {
        return CellRule::identifyFace<PyramidRule13>(element_nodes, nodes);
    }
}; /* class */

}  // namespace MeshLib
