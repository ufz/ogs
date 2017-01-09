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

class CellRule
{
public:
    /// Constant: Dimension of this mesh element
    static const unsigned dimension = 3u;

    /**
     * Checks if the node order of an element is correct by testing surface normals.
     * For 3D elements true is returned if the normals of all faces points away from the centre of
     * the element.
     * Note: This method might give wrong results if something else is wrong with the element
     * (non-planar faces, non-convex geometry, possibly zero volume) which causes the calculated
     * center of gravity to lie outside of the actual element
     */
    static bool testElementNodeOrder(const Element* /*e*/);
}; /* class */

} /* namespace */
