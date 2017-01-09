/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace MeshLib
{

class Element;

class VertexRule
{
public:
    /// Constant: Dimension of this mesh element
    static const unsigned dimension = 0u;

    /// Constant: The number of faces
    static const unsigned n_faces = 0;

    /// Constant: The number of edges
    static const unsigned n_edges = 0;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* /*e*/, unsigned /*i*/)
    {
        return nullptr;
    }

    /// Checks if the node order of an element is correct by testing surface
    /// normals.  For 0D elements this always returns true.
    static bool testElementNodeOrder(const Element* /*e*/) { return true; }

};

}
