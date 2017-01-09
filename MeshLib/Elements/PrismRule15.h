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
#include "PrismRule6.h"

namespace MeshLib
{

/**
 * This class represents a 3d prism element with 15 nodes. The following sketch shows the node and edge numbering.
 * @anchor PrismNodeAndEdgeNumbering
 * @code
 *            5
 *           / \
 *          / : \
 *       14/  :  \13
 *        /   :11 \
 *       /    : 12 \
 *      3-----------4
 *      |     :     |
 *      |     2     |
 *      |    . .    |
 *     9|   .   .   |10
 *      | 8.     .7 |
 *      | .       . |
 *      |.         .|
 *      0-----------1
 *            6
 *
 * @endcode
 */
class PrismRule15 : public PrismRule6
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 15u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::PRISM15;

    /// Constant: Local node index table for faces
    static const unsigned face_nodes[5][8];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[9][3];

    /// Constant: Table for the number of nodes for each face
    static const unsigned n_face_nodes[5];

    /// Returns the i-th edge of the element.
    typedef QuadraticEdgeReturn EdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

}; /* class */

} /* namespace */
