/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FACERULE_H_
#define FACERULE_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"

namespace MeshLib {

/**
 */
class FaceRule
{
public:
	/// Returns the face i of the element.
	static const Element* getFace(const Element* e, unsigned i) { return e->getEdge(i); }

	/// Get the number of nodes for face i.
	static unsigned getNFaceNodes(unsigned /*i*/) { return 2; }

	/// 2D elements have no faces.
	static unsigned getNFaces() { return 0; }
}; /* class */

} /* namespace */

#endif /* FACERULE_H_ */

