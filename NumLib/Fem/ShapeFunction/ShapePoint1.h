/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#include <array>

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

    static constexpr std::array reference_element_centre = {0.0};

    using MeshElement = MeshLib::Point;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
    static constexpr int ORDER = 0;
};
}  // namespace NumLib

#include "ShapePoint1-impl.h"
