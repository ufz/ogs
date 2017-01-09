/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FACERULE_H_
#define FACERULE_H_

#include "MathLib/Vector3.h"
#include "Element.h"

namespace MeshLib
{

class FaceRule
{
public:
    /// Constant: Dimension of this mesh element
    static const unsigned dimension = 2u;

    /// Returns the face i of the element.
    static const Element* getFace(const Element* e, unsigned i) { return e->getEdge(i); }

    /// Constant: The number of faces
    static const unsigned n_faces = 0;

    /**
     * Checks if the node order of an element is correct by testing surface normals.
     * For 2D elements true is returned if the normal points (roughly) upwards.
     */
    static bool testElementNodeOrder(const Element* /*e*/);

    /// Returns the surface normal of a 2D element.
    static MathLib::Vector3 getSurfaceNormal(const Element* e);

}; /* class */

} /* namespace */

#endif /* FACERULE_H_ */

