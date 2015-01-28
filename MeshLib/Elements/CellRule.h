/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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

}; /* class */

} /* namespace */

#endif /* HEXRULE_H_ */

