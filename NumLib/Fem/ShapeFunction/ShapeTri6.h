/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Tri.h"

namespace NumLib
{

/**
 *  Shape function for a 6-nodes triangle element
 *
 */
class ShapeTri6
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] N    a vector of calculated shape function.
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X &r, T_N &N);

    /**
     * Evaluate derivatives of the shape function at the given point
     * The point coordinates in \c r are not used.
     *
     * @param [in]  r    point coordinates
     * @param [out] dN   a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X &r, T_N &dN);

    using MeshElement = MeshLib::Tri6;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
};

}  // namespace NumLib

#include "ShapeTri6-impl.h"
