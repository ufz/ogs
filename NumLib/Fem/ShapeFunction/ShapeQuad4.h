/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once
#include <array>

#include "MeshLib/Elements/Quad.h"

namespace NumLib
{
/**
 *  Shape function for a quadrilateral element of four nodes in natural
 * coordinates
 *
 * \verbatim
 *  2 (-1,1)     1 (1,1)
 *     *--------*
 *     |        |
 *     |        |
 *     |        |
 *     *--------*
 *  3 (-1,-1)    4 (1,-1)
 * \endverbatim
 */
class ShapeQuad4
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] N   a vector of calculated shape function.
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X& r, T_N& N);

    /**
     * Evaluate derivatives of the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] dN  a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X& r, T_N& dN);

    static constexpr std::array reference_element_centre = {0.0, 0.0};

    using MeshElement = MeshLib::Quad;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
    static constexpr int ORDER = 1;
};

}  // namespace NumLib

#include "ShapeQuad4-impl.h"
