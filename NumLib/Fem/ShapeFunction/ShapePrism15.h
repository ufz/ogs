// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <array>

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

    static constexpr std::array reference_element_centre = {1. / 3., 1. / 3.,
                                                            0.0};

    using MeshElement = MeshLib::Prism15;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
    static constexpr int ORDER = 2;
};

}  // namespace NumLib

#include "ShapePrism15-impl.h"
