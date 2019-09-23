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

#include "MeshLib/Elements/Prism.h"

namespace NumLib
{

/**
 *  Shape function for a 15-nodes prism element in natural coordinates
 *
 */
class ShapePrism15
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  x   natural coordinates (r,s,t)
     * @param [out] N   a vector of calculated shape functions
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X& x, T_N& N);

    /**
     * Evaluate derivatives of the shape function at the given point
     *
     * @param [in]  x   natural coordinates (r,s,t)
     * @param [out] dN  a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X& x, T_N& dN);

    using MeshElement = MeshLib::Prism15;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
};

}  // namespace NumLib

#include "ShapePrism15-impl.h"
