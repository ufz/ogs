/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EDGERULE_H_
#define EDGERULE_H_

namespace MeshLib
{

class Element;

class EdgeRule
{
public:
    /// Constant: Dimension of this mesh element
    static const unsigned dimension = 1u;

    /// Constant: The number of faces
    static const unsigned n_faces = 0;

    /// Constant: The number of edges
    static const unsigned n_edges = 1;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* /*e*/, unsigned /*i*/) { return nullptr; }

    /**
    * Checks if the node order of an element is correct by testing surface normals.
    * For 1D elements this always returns true.
    */
    static bool testElementNodeOrder(const Element* /*e*/) { return true; }

}; /* class */

} /* namespace */

#endif /* EDGERULE_H_ */

