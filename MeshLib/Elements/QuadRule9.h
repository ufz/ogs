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
#include "QuadRule8.h"

namespace MeshLib
{

/**
 * This class represents a 2d quadrilateral element with 9 nodes.
 * The following sketch shows the node and edge numbering.
 * @anchor Quad9NodeAndEdgeNumbering
 * @code
 *              2
 *        3-----6-----2
 *        |           |
 *        |           |
 *      3 7     8     5 1
 *        |           |
 *        |           |
 *        0-----4-----1
 *              0
 * @endcode
 */
class QuadRule9 : public QuadRule8
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 9u;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::QUAD9;
}; /* class */

} /* namespace */
