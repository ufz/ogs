/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPEPRISM6_H_
#define SHAPEPRISM6_H_

#include "MeshLib/Elements/Prism.h"

namespace NumLib
{

/**
 *  Shape function for a 6-nodes prism element in natural coordinates
 *
 */
class ShapePrism6
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  r   natural coordinates (r,s,t)
     * @param [out] N   a vector of calculated shape functions
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X &r, T_N &N);

    /**
     * Evaluate derivatives of the shape function at the given point
     *
     * @param [in]  r   natural coordinates (r,s,t)
     * @param [out] dN  a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X &r, T_N &dN);

    using MeshElement = MeshLib::Prism;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
};

}

#include "ShapePrism6-impl.h"

#endif //SHAPEPRISM6_H_
