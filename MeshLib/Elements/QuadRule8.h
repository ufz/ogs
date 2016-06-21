/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUADRULE8_H_
#define QUADRULE8_H_

#include "MeshLib/MeshEnums.h"
#include "QuadRule4.h"
#include "EdgeReturn.h"

#ifdef _MSC_VER
#include "meshlib_export.h"
#else
#define MESHLIB_EXPORT
#endif

namespace MeshLib
{

/**
 * This class represents a 2d quadrilateral element with 8 nodes.
 * The following sketch shows the node and edge numbering.
 * @anchor Quad8NodeAndEdgeNumbering
 * @code
 *              2
 *        3-----6-----2
 *        |           |
 *        |           |
 *      3 7           5 1
 *        |           |
 *        |           |
 *        0-----4-----1
 *              0
 * @endcode
 */
class QuadRule8 : public QuadRule4
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 8u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::QUAD8;

    /// Constant: Local node index table for edge
    static MESHLIB_EXPORT const unsigned edge_nodes[4][3];

    /// Returns the i-th edge of the element.
    typedef QuadraticEdgeReturn EdgeReturn;

}; /* class */

} /* namespace */

#endif /* QUADRULE8_H_ */

