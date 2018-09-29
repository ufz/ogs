/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "EdgeReturn.h"
#include "TetRule4.h"

namespace MeshLib
{

/**
 * This class represents a 3d tetrahedron element with 10 nodes. The following sketch shows the node and edge numbering.
 * @anchor TetrahedronNodeAndEdgeNumbering
 * @code
 *          3
 *         /|\
 *        / | \
 *      7/  |  \9
 *      /   |8  \
 *     /    |    \
 *    0.....|.....2
 *     \    |  2 /
 *      \   |   /
 *      4\  |  /5
 *        \ | /
 *         \|/
 *          1
 *
 * @endcode
 */
class TetRule10 : public TetRule4
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 10u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::TET10;

    /// Constant: Local node index table for faces
    static const unsigned face_nodes[4][6];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[6][3];

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

}; /* class */

} /* namespace */
