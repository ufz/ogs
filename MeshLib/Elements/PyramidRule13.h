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
#include "Element.h"
#include "EdgeReturn.h"
#include "PyramidRule5.h"

namespace MeshLib
{

/**
 * This class represents a 3d pyramid element. The following sketch shows the node and edge numbering.
 * @anchor PyramidNodeAndEdgeNumbering
 * @code
 *
 *               4
 *             //|\
 *            // | \
 *         12//  |  \11
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
 * @endcode

 */
class PyramidRule13 : public PyramidRule5
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
    typedef QuadraticEdgeReturn EdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

}; /* class */

} /* namespace */
