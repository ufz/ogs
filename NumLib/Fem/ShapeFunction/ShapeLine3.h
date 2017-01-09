/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPELINE3_H_
#define SHAPELINE3_H_

#include "MeshLib/Elements/Line.h"

namespace NumLib
{

/**
 *  Shape function for a 3-nodes line element
 */
class ShapeLine3
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] N   a vector of calculated shape function.
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X &r, T_N &N);

    /**
     * Evaluate derivatives of the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] dN  a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X &r, T_N &dN);

    using MeshElement = MeshLib::Line3;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
};

}

#include "ShapeLine3-impl.h"

#endif //SHAPELINE3_H_
