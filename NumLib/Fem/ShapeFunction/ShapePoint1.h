/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPEPOINT1_H_
#define SHAPEPOINT1_H_

#include "MeshLib/Elements/Point.h"

namespace NumLib
{
///  Shape function for a point element in natural coordinates.
class ShapePoint1
{
public:
    /// Evaluate the shape function at the given point
    ///
    /// @param [in]  r    point coordinates
    /// @param [out] N   a vector of calculated shape function.
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X& r, T_N& N);

    using MeshElement = MeshLib::Point;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
};
}

#include "ShapePoint1-impl.h"

#endif  // SHAPEPOINT1_H_
