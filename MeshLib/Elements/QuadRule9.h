/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUADRULE9_H_
#define QUADRULE9_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "QuadRule8.h"

namespace MeshLib
{

/**
 * This class represents a 2d quadliteral element. The following sketch shows the node and edge numbering.
 * @anchor QuadNodeAndEdgeNumbering
 * @code
 *              2
 *        3-----------2
 *        |           |
 *        |           |
 *       3|           |1
 *        |           |
 *        |           |
 *        0-----------1
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

#endif /* QUADRULE9_H_ */

