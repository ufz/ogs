/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the TemplateHex class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CELLRULE_H_
#define CELLRULE_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"

namespace MeshLib {

/**
 */
class CellRule
{
public:
	/// Constant: Dimension of this mesh element
	static const unsigned dimension = 3u;

	/**
	 * Checks if the node order of an element is correct by testing surface normals.
	 */
	static bool testElementNodeOrder(const Element* e);

}; /* class */

} /* namespace */

#endif /* HEXRULE_H_ */

