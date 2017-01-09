/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef HEXRULE20_H_
#define HEXRULE20_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "EdgeReturn.h"
#include "HexRule8.h"

namespace MeshLib
{

/**
 * A 20-nodes Hexahedron Element.
 * @code
 *
 *  Hex:
 *                14
 *          7-----------6
 *         /:          /|
 *        / :         / |
 *     15/  :        /13|
 *      / 19:       /   | 18
 *     /    : 12   /    |
 *    4-----------5     |
 *    |     :     | 10  |
 *    |     3.....|.....2
 *    |    .      |    /
 * 16 |   .       |17 /
 *    |11.        |  / 9
 *    | .         | /
 *    |.          |/
 *    0-----------1
 *          8
 *
 * @endcode
 */
class HexRule20 : public HexRule8
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 20u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::HEX20;

    /// Constant: Local node index table for faces
    static const unsigned face_nodes[6][8];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[12][3];

    /// Returns the i-th edge of the element.
    typedef QuadraticEdgeReturn EdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

}; /* class */

} /* namespace */

#endif /* HEXRULE20_H_ */

